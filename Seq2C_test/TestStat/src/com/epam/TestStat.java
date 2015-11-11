
package com.epam;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Scanner;
import org.apache.commons.math.stat.descriptive.moment.Mean;
import org.apache.commons.math.stat.descriptive.rank.Max;
import org.apache.commons.math.stat.descriptive.rank.Median;
import org.apache.commons.math.stat.descriptive.rank.Min;
import org.apache.commons.math.stat.descriptive.rank.Percentile;
import org.apache.commons.math.stat.descriptive.summary.Sum;

/**
 *
 * @author Sofya Zhuk
 */
public class TestStat {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws FileNotFoundException {
        Date dStart = new Date();
        List<Double> list = new ArrayList<>();
        File file=new File(args[0]);
        Scanner in = new Scanner(file);
        while (in.hasNext()){
            double d = Double.parseDouble(in.next());
            list.add(d);
        }
        double[] initialDoubleArray=new double[list.size()];    
        for(int i=0;i<=list.size()-1;i++){
            initialDoubleArray[i]=list.get(i);
        }
        double result;
        Date dEnd = new Date();
        System.out.println("Time for preparation: "+(dEnd.getTime()-dStart.getTime())+" ms");
        //test for mean
        Mean mean=new Mean();
        Date d0 = new Date();
        result=mean.evaluate(initialDoubleArray);
        Date d1 = new Date();
        System.out.println(result);
        System.out.println("Time for mean: "+(d1.getTime()-d0.getTime())+" ms");
        //test for median
        Median median=new Median();
        Date d2 = new Date();
        result=median.evaluate(initialDoubleArray);
        Date d3 = new Date();
        System.out.println(result);
        System.out.println("Time for median: "+(d3.getTime()-d2.getTime())+" ms");
        //test for max
        Max max=new Max();
        Date d4 = new Date();
        result=max.evaluate(initialDoubleArray);
        Date d5 = new Date();
        System.out.println(result);
        System.out.println("Time for max: "+(d5.getTime()-d4.getTime())+" ms");
        //test for min
        Min min=new Min();
        Date d6 = new Date();
        result=min.evaluate(initialDoubleArray);
        Date d7 = new Date();
        System.out.println(result);
        System.out.println("Time for min: "+(d7.getTime()-d6.getTime())+" ms");
        //test for sum
        Sum sum=new Sum();
        Date d8 = new Date();
        result=sum.evaluate(initialDoubleArray);
        Date d9 = new Date();
        System.out.println(result);
        System.out.println("Time for sum: "+(d9.getTime()-d8.getTime())+" ms");
        //test for percentile
        Percentile percentile=new Percentile();
        Date d10 = new Date();
        result=percentile.evaluate(initialDoubleArray, 10);
        Date d11 = new Date();
        System.out.println(result); 
        System.out.println("Time for percentile: "+(d11.getTime()-d10.getTime())+" ms");   
    }
    
}
