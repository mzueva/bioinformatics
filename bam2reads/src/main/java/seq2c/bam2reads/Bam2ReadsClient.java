package seq2c.bam2reads;


/**
 * Client for launching Bam2Reads tool
 * Reads two arguments from command line
 * fist argument = file, containing list of bam files
 * second argument = output file for results
 * example of usage: bam2reads.jar sample2bam.txt stat.txt
 */

public class Bam2ReadsClient {
    public static void main(String[] args) {
        if (args.length == 2) {
            Bam2Reads.printStatsToFile(args[0], args[1]);
        } else {
            System.err.println("Bam2Reads needs two arguments: input file and output file");
        }
    }
}
