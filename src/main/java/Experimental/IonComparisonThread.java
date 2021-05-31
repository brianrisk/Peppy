package Experimental;

import Peppy.Peptide;
import Peppy.Properties;

import java.util.ArrayList;
import java.util.Hashtable;

public class IonComparisonThread implements Runnable {

    ArrayList<Peptide> peptides;
    IonComparisonServer ionComparisonServer;

    /**
     * @param proteins
     * @param spectrum
     */
    public IonComparisonThread(ArrayList<Peptide> peptides, IonComparisonServer ionComparisonServer) {
        this.peptides = peptides;
        this.ionComparisonServer = ionComparisonServer;

    }


    public void run() {
        while (peptides != null) {
            Peptide peptideA, peptideB;
            for (int indexOne = 0; indexOne < peptides.size() - 1; indexOne++) {
                peptideA = peptides.get(indexOne);
                int peptideMassIndex = (int) (Math.round(peptideA.getMass()) - Properties.peptideMassMinimum);
                Hashtable<Integer, Integer> peptideAIons = CommonIonCount.getIons(peptideA);
                for (int indexTwo = indexOne + 1; indexTwo < peptides.size(); indexTwo++) {
                    peptideB = peptides.get(indexTwo);
                    if (peptideA.equals(peptideB)) continue;
                    Hashtable<Integer, Integer> peptideBIons = CommonIonCount.getIons(peptideB);
                    int numberOfPeaksInCommon = CommonIonCount.compare(peptideMassIndex, peptideAIons, peptideBIons);
                }
            }

            peptides = ionComparisonServer.getNextObject();
        }
    }


}