package Experimental;

import Peppy.U;

import java.util.Hashtable;
import java.util.Random;

public class JoshProbabilityExperiment {

    public static void main(String args[]) {
        Random random = new Random();
        Hashtable<Integer, Integer> totals = new Hashtable<Integer, Integer>();

        for (int experimentCount = 0; experimentCount < 100000; experimentCount++) {

            Hashtable<Integer, Boolean> usedNumbers = new Hashtable<Integer, Boolean>();

            int max = 100;

            int repeatCount = 0;
            int uniqueCount = 0;
            while (uniqueCount < 100) {
                int grabbed = random.nextInt(max);
                Boolean previouslyGrabbed = usedNumbers.get(grabbed);
                if (previouslyGrabbed == null) {
                    usedNumbers.put(grabbed, true);
                    uniqueCount++;
                } else {
                    repeatCount++;
                }
            }

            /* now that we have the run count, add it to the toal */
            Integer total = totals.get(repeatCount);
            if (total == null) {
                totals.put(repeatCount, 1);
            } else {
                totals.put(repeatCount, total + 1);
            }
        }

        /* print all the contents of totals */
        for (int key : totals.keySet()) {
            int total = totals.get(key);
            U.p(key + "\t" + total);
        }

    }

}
