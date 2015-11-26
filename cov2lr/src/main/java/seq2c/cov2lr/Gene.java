package seq2c.cov2lr;

/**
 * Created by Mariia_Zueva on 11/25/2015.
 */
public class Gene {
    private String name;
    private String chr;
    private long start;
    private long end;
    private long len;

    public Gene(String name, String chr, long start, long end, long len) {
        this.name = name;
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.len = len;
    }

    public long getStart() {
        return start;
    }

    public void setStart(long start) {
        this.start = start;
    }

    public long getEnd() {
        return end;
    }

    public void setEnd(long end) {
        this.end = end;
    }

    public long getLen() {
        return len;
    }

    public String getName() {
        return name;
    }

    public String getChr() {
        return chr;
    }

    public void setLen(long len) {
        this.len = len;
    }

    @Override
    public String toString() {
        return name + " " + chr + " " + start + " " + end + " " + len;
    }
}

