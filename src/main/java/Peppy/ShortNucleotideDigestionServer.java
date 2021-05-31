package Peppy;

import java.util.ArrayList;

/**
 * Copyright 2013, Brian Risk
 *
 * @author Brian Risk
 */
public class ShortNucleotideDigestionServer {

    ArrayList<NucleotideSequence> nucleotideSequences;

    /* every nucleotideSequence produces a "cluster" of proteins */
    ArrayList<ArrayList<Protein>> proteinClusters;

    ArrayList<Thread> threads;

    /* this is how we keep track of which protein to give out next */
    private int nucleotideIndex = 0;

    /* why this simply isn't Properties.numberOfThreads is because we may have less nucleotideSequences than that number */
    private int numberOfThreads;


    public ShortNucleotideDigestionServer(ArrayList<NucleotideSequence> nucleotideSequences, boolean reverseDatabase) {
        this.nucleotideSequences = nucleotideSequences;
        proteinClusters = new ArrayList<ArrayList<Protein>>(nucleotideSequences.size());

        //here we make sure we don't use more threads than we have proteins
        numberOfThreads = Properties.numberOfThreads;
        if (Properties.numberOfThreads > nucleotideSequences.size()) numberOfThreads = nucleotideSequences.size();
        threads = new ArrayList<Thread>(numberOfThreads);

        //spawn new threads as needed
        for (int threadNumber = 0; threadNumber < numberOfThreads; threadNumber++) {
            ShortNucleotideDigestionThread digester = new ShortNucleotideDigestionThread(getNextElementToProcess(), this, reverseDatabase);
            Thread thread = new Thread(digester);
            thread.start();
            threads.add(thread);
        }
    }


    public synchronized void report(ArrayList<Protein> proteinsFromOneNucleotideSequence) {
        proteinClusters.add(proteinsFromOneNucleotideSequence);
    }


    public synchronized NucleotideSequence getNextElementToProcess() {
        NucleotideSequence out = null;
        if (nucleotideIndex < nucleotideSequences.size()) {
            out = nucleotideSequences.get(nucleotideIndex);
            nucleotideIndex++;
        }
        return out;
    }

    /**
     * First we wait for all of  our  threads to finish.
     *
     * @return accumulated
     */
    public ArrayList<Protein> getResults() {
        boolean going = true;
        while (going) {
            for (int threadNumber = 0; threadNumber < numberOfThreads; threadNumber++) {
                Thread thread = threads.get(threadNumber);
                going = thread.isAlive();
                //if at least one thread is going, break out of this for loop
                if (going) break;
            }
            //Sleep for a bit and wait for the threads to finish.
            try {
                Thread.sleep(500);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
        //calculate size of combined ArrayLists
        int size = 0;
        for (ArrayList<Protein> proteinCluster : proteinClusters) {
            size += proteinCluster.size();
        }
        //combine matches into result and return
        ArrayList<Protein> out = new ArrayList<Protein>(size);
        for (ArrayList<Protein> proteinCluster : proteinClusters) {
            out.addAll(proteinCluster);
        }
        return out;
    }

}
