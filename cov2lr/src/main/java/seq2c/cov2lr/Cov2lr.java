package seq2c.cov2lr;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 *  Normalize the coverage from targeted sequencing to CNV log2 ratio.  The algorithm assumes the medium
 *     is diploid, thus not suitable for homogeneous samples (e.g. parent-child).
 */
public class Cov2lr {

    /**
      *  Indicate this is amplicon or exon based calling.  By default, it'll aggregate at gene level.
     */
    private boolean amplicon;

    /**
     * Use the control sample(s)
     */
    private boolean useControlSamples;

    /**
      *  Array of control samples
     */
    private String[] controlSamples;


    /**
     *   Statistics from first input file
     *   Key = name of sample
     *   Value = number of aligned reads
     */
    private Map<String, Long> mappingReads;

    /**
      *  Normalization factor for each sample
      *  calculated as mean number of reads / number of reads per sample
      *  Key = name of sample
      *  Value = factor
     */
    private Map<String, Double> factor;
    /**
     *  Map of genes from second input file
     *  Key = name of gene
     *  Value = gene object
     */
    private Map<String, Gene> genes;
    /**
     *  Map of samples
     *  Key = key for string, consists of sample name, gene name, start, end and length
     *  Value = map of sample names and sample objects
     */
    private Map<String, Map<String, Sample>> samples;

    /**
     *  The failed factor for individual amplicons.
     *  If (the 80th percentile of an amplicon depth)/(the global median depth)
     *  is less than the argument, the amplicon is considered failed and won't be used in calculation.  Default: 0.2.
     */
    private final double FAILEDFACTOR = 0.2;

    /**
     * Constructor reads the input files and constructs genes, samples, factor maps
     * @param amplicon  =   determines the aggregation level (gene or record)
     * @param fileStat  =   A file containing # of mapped or sequenced reads for samples.  At least two columns.
     *                      First is the sample name, 2nd is the number of mapped or sequenced reads.
     * @param fileCov   =   The coverage output file from checkCov.pl script.  Can also take from standard in or more than
     *                      one file.
     * @param useControlSamples = use the control sample(s)
     * @param controlSamples     = multiple controls are allowed, which are separated by ":"
     *
     */

    public Cov2lr(boolean amplicon, String fileStat, String fileCov, boolean useControlSamples, String controlSamples) {
        this.amplicon = amplicon;
        this.useControlSamples = false;
        this.mappingReads = readStat(fileStat);
        this.genes = new LinkedHashMap<>();
        this.samples = new LinkedHashMap<>();
        readCoverage(fileCov, this.genes, this.samples);
        this.factor = new LinkedHashMap<>();
        setFactor();
        this.useControlSamples = useControlSamples;
        if (useControlSamples) {
            this.controlSamples = controlSamples.split(":");
        }
    }


    /**
     * Main method, makes the statistics calculation according to the algorithm
     * print result to the Standart output
     */

    public void doWork() {
        double[] depth = new double[samples.size()];

        Set<String> samp = new HashSet<>();

        double medDepth = getMedDepth(depth, samp);

        List<String> bad = new LinkedList<>();
        List<Double> gooddepth = new LinkedList<>();

        splitQualitySamples(samp, medDepth, bad, gooddepth);

        medDepth = new Median().evaluate(toDoubleArray(gooddepth));

        Map<String, Double> factor2 = new LinkedHashMap<>();

        setFactor2(medDepth, bad, factor2);

        Map<String, Double> sampleMedian = getSampleMedian(samp, bad);

        setNorm(medDepth, bad, factor2, sampleMedian);

        printResult(bad);
    }


    private Map<String, Long> readStat(String fileName) {
        Map<String, Long> map = new LinkedHashMap<>();
        try (BufferedReader reader= new BufferedReader(new FileReader(fileName))) {
            String line;
            while ((line = reader.readLine()) != null) {
                String[] samples = line.split("\\s+");
                if (samples.length == 2) {
                    map.put(samples[0], Long.parseLong(samples[1]));
                }
            }
        } catch (IOException e) {
            System.err.println("Cannot open file " + fileName);
        } catch (NumberFormatException e) {
            System.err.println("Cannot parse long number");
        }
        return map;
    }

    private void readCoverage(String fileName, Map<String, Gene> genes, Map<String, Map<String, Sample>> samples) {
        try (BufferedReader reader= new BufferedReader(new FileReader(fileName))) {
            String line;
            while ((line = reader.readLine()) != null) {
                //skip lines with words "Sample, Whole, Control, Undertermined
                if (line.contains("Sample") || line.contains("Whole") || line.contains("Control") ||
                        line.contains("Undertermined")) continue;

                String[] sampleLines = line.split("\\s+");
                if (sampleLines.length == 8) {

                    String sample = sampleLines[0];
                    String gene = sampleLines[1];
                    String chr = sampleLines[2];
                    long start = Long.parseLong(sampleLines[3]);
                    long end = Long.parseLong(sampleLines[4]);
                    long len = Long.parseLong(sampleLines[6]);
                    long depth = Long.parseLong(sampleLines[7]);

                    addGene(genes, gene, chr, start, end, len);
                    String key = amplicon ? concatKey(gene, chr, sampleLines[3], sampleLines[4], sampleLines[6]) : gene ;
                    addSample(samples, key, sample, len, depth, gene);
                }
            }
        } catch (IOException e) {
            System.err.println("Cannot open file " + fileName);
        } catch (NumberFormatException e) {
            System.err.println("Cannot parse long number");
        }

    }

    private void addGene(Map<String, Gene> map, String gene, String chr, long start, long end, long len) {


        if (map.containsKey(gene)) {
            Gene geneObj = map.get(gene);
            //move start
            if (geneObj.getStart() > start) {
                geneObj.setStart(start);
            }
            //move end
            if (geneObj.getEnd() < end) {
                geneObj.setEnd(end);
            }
            //update length
            geneObj.setLen(geneObj.getLen() + len);
        } else {
            Gene geneObj = new Gene(gene, chr, start, end, len);
            map.put(gene, geneObj);
        }
    }

    private void addSample(Map<String, Map<String, Sample>> map, String key, String sample, long len, long depth, String gene) {

        if (map.containsKey(key)) {
            Map<String, Sample> sampleMap = map.get(key);
            if (sampleMap.containsKey(sample)) {
                Sample sampleObj = sampleMap.get(sample);
                sampleObj.addLen(len);
                sampleObj.addCov(depth);
            } else {
                Sample sampleObj = new Sample(key, sample, gene, len, depth);
                sampleMap.put(sample, sampleObj);
            }
        } else {
            Sample sampleObj = new Sample(key, sample, gene, len, depth);
            Map<String, Sample> sampleMap = new LinkedHashMap<>();
            sampleMap.put(sample, sampleObj);
            map.put(key, sampleMap);

        }
    }

    private void setFactor() {
        double[] array = new double[mappingReads.size()];
        int i = 0;
        for (Long l : mappingReads.values()) {
            array[i++] = l;
        }
        double meanReads = new Mean().evaluate(array);

        for (Map.Entry<String, Long> entry : mappingReads.entrySet()) {
            factor.put(entry.getKey(), meanReads / entry.getValue());
        }

    }

    public void doTest() {

        try (FileWriter writer = new FileWriter("java_debug.txt")) {
            writer.write("Factor:" + "\n");
            for (Map.Entry<String, Double> entry: factor.entrySet()) {
                writer.write(entry.getKey() + "\t" + entry.getValue() + "\n");
                writer.flush();
            }
            writer.write("*********************" + "\n");
            writer.write("Loc hash:" + "\n");
            writer.flush();
            long start = 0;
            long end = 0;
            long len = 0;
            for (Map.Entry<String, Gene> entry : genes.entrySet()) {
                start += entry.getValue().getStart();
                end += entry.getValue().getEnd();
                len += entry.getValue().getLen();

                //writer.write(entry.getKey() + "\t" + entry.getValue() + "\n");
                //writer.flush();
            }

            writer.write("start = " + start + " ");
            writer.write("end = " + end + " ");
            writer.write("len = " + len + "\n");
            writer.write("*********************" + "\n");
            writer.write("Data hash:" + "\n");
            writer.flush();
            long n = samples.size();
            long lenS = 0;
            long cov = 0;
            for (Map.Entry<String, Map<String, Sample>> entry : samples.entrySet()) {
                for (Map.Entry<String, Sample> sampleEntry : entry.getValue().entrySet()) {
                    lenS += sampleEntry.getValue().getLen();
                    cov += sampleEntry.getValue().getCov();
                }
            }
            writer.write("data len = " + lenS + " ");
            writer.write("cov = " + cov + " ");
            writer.write("numK = " + n + "\n");
            writer.write("*********************" + "\n");

        } catch (IOException e) {
            System.err.println(e.getLocalizedMessage());
        }
    }

    private void printResult(List<String> bad) {
        for (Map.Entry<String, Map<String, Sample>> entry : samples.entrySet()) {
            if (!bad.contains(entry.getKey())) {
                for (Map.Entry<String, Sample> entrySample : entry.getValue().entrySet()) {
                    Sample sample = entrySample.getValue();
                    String title = amplicon ? entry.getKey() : sample.getTitle(genes.get(sample.getGene()));
                    String result = sample.getResultString(title);
                    if (useControlSamples) {
                        addControlSamples(result, sample);
                    }
                    System.out.println(result);
                }
            }
        }
    }

    private void addControlSamples(String result, Sample sample) {
        List<Double> list = new LinkedList<>();
        for (String s: controlSamples) {
            for (Map.Entry<String, Map<String, Sample>> entry : samples.entrySet()) {
                list.add(entry.getValue().get(s).getNorm1b());
            }
        }
        if (!list.isEmpty()) {
            double mean = new Mean().evaluate(toDoubleArray(list));
            StringBuilder builder = new StringBuilder(result);
            builder.append("\t").append(String.format("%.3f%n", mean == 0 ? sample.getNorm1b() / mean / Math.log(2) : 0));
        }
    }

    private void setNorm(double medDepth, List<String> bad, Map<String, Double> factor2, Map<String, Double> sampleMedian) {
        for (Map.Entry<String, Map<String, Sample>> entry : samples.entrySet()) {
            if (!bad.contains(entry.getKey())) {
                for (Map.Entry<String, Sample> entrySample : entry.getValue().entrySet()) {
                    Sample sample = entrySample.getValue();
                    double norm1 = sample.getNorm1();
                    double fact2 = factor2.get(entry.getKey());
                    double smplMed = sampleMedian.get(sample.getSample());
                    sample.setNorm1b(norm1 * fact2 + 0.1);
                    sample.setNorm2(round(medDepth != 0 ? Math.log((norm1 * fact2 + 0.1) / medDepth / Math.log(2)) : 0));
                    sample.setNorm3(round(smplMed != 0 ? Math.log((norm1 * fact2 + 0.1) / smplMed / Math.log(2)) : 0));

                }
            }
        }
    }

    private Map<String, Double> getSampleMedian(Set<String> samp, List<String> bad) {
        Map<String, Double> sampleMedian = new LinkedHashMap<>();

        for (String s : samp) {
            List<Double> list = new LinkedList<>();
            for (Map.Entry<String, Map<String, Sample>> entry : samples.entrySet()) {
                if (!bad.contains(entry.getKey())) {
                    list.add(entry.getValue().get(s).getNorm1());
                }
            }
            sampleMedian.put(s, new Median().evaluate(toDoubleArray(list)));
        }
        return sampleMedian;
    }

    private void setFactor2(double medDepth, List<String> bad, Map<String, Double> factor2) {
        for (Map.Entry<String, Map<String, Sample>> entry : samples.entrySet()) {
            if (!bad.contains(entry.getKey())) {
                List<Double> list = new LinkedList<>();
                for (Map.Entry<String, Sample> entrySample : entry.getValue().entrySet()) {
                    Sample sample = entrySample.getValue();
                    list.add(sample.getNorm1());
                }
                double median = new Median().evaluate(toDoubleArray(list));
                if (median != 0) {
                    factor2.put(entry.getKey(), medDepth/median);
                } else {
                    factor2.put(entry.getKey(), 0.0);
                }
            }
        }
    }

    private void splitQualitySamples(Set<String> samp, double medDepth, List<String> bad, List<Double> gooddepth) {
        for (Map.Entry<String, Map<String, Sample>> entry : samples.entrySet()) {
            List<Double> temp = new LinkedList<>();
            double kp80 = filterData(samp, entry.getValue(), temp);
            if (kp80 < medDepth * FAILEDFACTOR) {
                bad.add(entry.getKey());
            } else {
                gooddepth.addAll(temp);
            }
        }
    }

    private double getMedDepth(double[] depth, Set<String> samp) {
        int i = 0;
        for (Map.Entry<String, Map<String, Sample>> entry: samples.entrySet()) {
            Map<String, Sample> sampleMap =  entry.getValue();

            for(Map.Entry<String, Sample> sampleEntry : sampleMap.entrySet()) {
                Sample sample = sampleEntry.getValue();
                double norm1 = round(sample.getCov() * factor.get(sample.getSample()));
                sample.setNorm1(norm1);
                depth[i++] = norm1;
                samp.add(sample.getSample());
            }
        }
        return new Median().evaluate(depth);
    }

    private double round(double d) {
        return Math.round(d * 100) / 100;
    }

    private double filterData(Set<String> sampleSet, Map<String, Sample> map, List<Double> list) {

        for (String sample : sampleSet) {
            if (map.containsKey(sample)) {
                list.add(map.get(sample).getNorm1());
            }
        }
        double[] result = toDoubleArray(list) ;
        return new Percentile().evaluate(result, 80);
    }

    private double[] toDoubleArray(List<Double> list) {
        double[] array = new double[list.size()];
        int i = 0;
        for (Double d : list) {
            array[i++] = d;
        }
        return array;
    }

    private String concatKey(String ... array) {
        StringBuilder builder = new StringBuilder();
        for (String s : array) {
            builder.append(s).append(" ");
        }
        return builder.toString().trim();
    }

}
