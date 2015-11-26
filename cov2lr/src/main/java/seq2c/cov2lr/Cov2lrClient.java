package seq2c.cov2lr;

/**
 * Created by Mariia_Zueva on 11/25/2015.
 */
public class Cov2lrClient {
    public static void main(String[] args) {
        if (args.length == 2) {
            boolean a = true; //amplicon or exon based calling
            boolean c = true; //use samples
            String samples = "AURA_26";
            Cov2lr script = new Cov2lr(a, args[0], args[1], c, samples);
            script.doWork();
        }
    }
}
