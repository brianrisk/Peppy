package Experimental;

import Peppy.Peptide;
import Peppy.Properties;
import Peppy.U;

import java.util.ArrayList;

public class IonComparisonServer {


    ArrayList<ArrayList<Peptide>> collectedArrays;

    //this is how we keep track of which Spectrum to give out next
    private int index = 0;

    ArrayList<Thread> threads;

    /* why this simply isn't Properties.numberOfThreads is because we may have less spectra than that number */
    private int numberOfThreads;


    public IonComparisonServer(ArrayList<ArrayList<Peptide>> collectedArrays) {
        this.collectedArrays = collectedArrays;


        //here we make sure we don't use more threads than we have proteins
        numberOfThreads = Properties.numberOfThreads;
        if (Properties.numberOfThreads > collectedArrays.size()) numberOfThreads = collectedArrays.size();
        threads = new ArrayList<Thread>(numberOfThreads);

        //spawn new threads as needed
        for (int threadNumber = 0; threadNumber < numberOfThreads; threadNumber++) {
            IonComparisonThread threadObject = new IonComparisonThread(getNextObject(), this);
            Thread thread = new Thread(threadObject);
            thread.start();
            threads.add(thread);
        }
    }


    /**
     * takes no parameters because it is used at the start when there are no matches to incorporate.
     *
     * @return
     */
    public synchronized ArrayList<Peptide> getNextObject() {
        ArrayList<Peptide> out = null;
        if (index < collectedArrays.size()) {
            out = collectedArrays.get(index);
            index++;
        }
        return out;
    }

    /**
     * First we wait for all of  our threads to finish.
     * "Yo threads, I'm happy for you and I'ma let you finish"
     *
     * @return accumulated
     */
    public void findMatches() {
        boolean going = true;
        while (going) {
            for (Thread thread : threads) {
                going = thread.isAlive();
                //if at least one thread is going, break out of this for loop
                if (going) break;
            }
            //Sleep for a bit and wait for the threads to finish.
            try {
                Thread.sleep(500);
            } catch (InterruptedException e) {
                U.p("getMatches in ScoringEngine thread was interrupted!");
                e.printStackTrace();
            }
        }

    }

}
