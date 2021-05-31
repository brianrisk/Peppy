package Peppy;

import java.util.ArrayList;

/**
 * Copyright 2013, Brian Risk
 *
 * @author Brian Risk
 */
public class VariableModificationServer {

    ArrayList<Peptide> peptides;

    ArrayList<Thread> threads;


    /**
     * takes no parameters because it is used at the start when there are no matches to incorporate.
     *
     * @return
     */
    public synchronized MatchesSpectrum getNextPeptide() {
        MatchesSpectrum out = null;
//		if (spectrumMatchesIndex < spectraMatches.size()) {
//			out = spectraMatches.get(spectrumMatchesIndex);
//			spectrumMatchesIndex++;
//		}
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
