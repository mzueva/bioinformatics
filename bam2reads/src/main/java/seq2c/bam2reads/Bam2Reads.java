package seq2c.bam2reads;

import htsjdk.samtools.AbstractBAMFileIndex;
import htsjdk.samtools.BAMIndexMetaData;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.*;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Tool for counting number of aligned records for list of bam files, given in a text file.
 * Tool reads the number of aligned records for each bam file from its index
 * and writes result into an output file.
 */
public class Bam2Reads {

    /**
     * Reads a map of sample name/bam filename from fileIn. For each bam file reads its indices,
     * using samtools and gets a number of aligned records from index information.
     * Number of aligned records for each file is recorded into output file,
     * Format of output file:
     * one line per bam file
     * each line: "sample name" \t "number of aligned records"
     *
     * @param fileIn  contains list of sample names and bam files to process,
     *                each line in file represent one file
     * @param fileOut output for results
     */

    public static void printStatsToFile(String fileIn, String fileOut) {
        Map<String, String> files = parseFile(fileIn);
        try (FileWriter writer = new FileWriter(fileOut)) {
            for (Map.Entry<String, String> entry : files.entrySet()) {
                writer.write(entry.getKey() + "\t");
                try (SamReader sam = SamReaderFactory.makeDefault().open(new File(entry.getValue()))) {

                    SamReader.Indexing ind = sam.indexing();
                    AbstractBAMFileIndex index = (AbstractBAMFileIndex) ind.getIndex();
                    int count = 0;
                    for (int i = 0; i < index.getNumberOfReferences(); i++) {
                        BAMIndexMetaData meta = index.getMetaData(i);
                        count += meta.getAlignedRecordCount();
                    }
                    writer.write(count + "\n");
                    writer.flush();
                } catch (IOException e) {
                    System.err.println("Cannot read file " + entry.getValue());
                }
            }
        } catch (IOException e) {
            System.err.println("Cannot write to file " + fileOut);
        }
    }

    /**
     * Reads input file line by line and creates a map containing
     * file header/bam file name
     *
     * @param fileName contains list of sample names and bam files to process,
     *                 each line in file represent one file
     * @return map with bam file names
     */
    private static Map<String, String> parseFile(String fileName) {
        Map<String, String> map = new LinkedHashMap<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(new File(fileName)))) {
            String currentLine;
            while ((currentLine = reader.readLine()) != null) {
                String[] samples = currentLine.split("\\s+");
                if (samples.length == 2) {
                    map.put(samples[0], samples[1]);
                }
            }
        } catch (IOException e) {
            System.err.println("Cannot open file " + fileName);
        }
        return map;
    }

}


