package Peppy;

import java.util.ArrayList;

/**
 * A thread that is given a peptide, and a list of variable modifications and, in return
 * produces a list of all potential modifications of that peptide
 * <p>
 * Copyright 2013, Brian Risk
 *
 * @author Brian Risk
 */
public class VariableModificationThread implements Runnable {

    /* the peptide to which we will find modifications */
    Peptide peptide;

    /* our list of potential modifications */
    ArrayList<ModificationVariable> modificaitons;

    /* our commanding server; results get reported back to this and new spectra to search come from here */
    VariableModificationServer variableModificationServer;

    public VariableModificationThread(Peptide peptide,
                                      ArrayList<ModificationVariable> modificaitons,
                                      VariableModificationServer variableModificationServer) {
        super();
        this.peptide = peptide;
        this.modificaitons = modificaitons;
        this.variableModificationServer = variableModificationServer;
    }


    public void run() {

        while (peptide != null) {

            ArrayList<PeptideWithModifications> modifiedPeptideList = new ArrayList<PeptideWithModifications>();
            byte[] sequence = peptide.getAcidSequence();
            String peptideString = peptide.getAcidSequenceString();

            int[] modificationIndices = new int[sequence.length];
            for (int index = 0; index < modificationIndices.length; index++) {
                modificationIndices[index] = -1;
            }

            /*
             * This function returns the acid being modified, but also updates
             * modificationIndices to correspond to the next viable modification
             */
            int residueIndex = getNextRadixArray(sequence, 0, modificationIndices);

            while (residueIndex != -1) {
                double[] modificaitonArray = new double[sequence.length];
                for (int index = 0; index < modificaitonArray.length; index++) {
                    if (modificationIndices[index] == -1) {
                        modificaitonArray[index] = 0;
                    } else {
                        if (sequence[index] == modificaitons.get(modificationIndices[index]).getAminoAcid()) {
                            modificaitonArray[index] = modificaitons.get(modificationIndices[index]).getMass();
                        }
                    }
                }
                modifiedPeptideList.add(new PeptideWithModifications(peptideString, modificaitonArray));


                residueIndex = getNextRadixArray(sequence, 0, modificationIndices);
            }

            /* return results, get new task */
//			peptide = variableModificationServer.getNextPeptide(modifiedPeptideList);

        }
    }


    private int getNextRadixArray(byte[] sequence, int residueIndex, int[] modificationIndices) {

        /* keep going until we find mod configuration that actually applies to the amino acids*/
        boolean keepGoing = true;

        /* if we went from one digit to the next in this step */
        boolean carryTheZero = false;

        while (keepGoing) {
            /* increment the indicies being sure to carry the zero */
            modificationIndices[residueIndex]++;

            /* e.g. 9999 to 10000 loop would repeat 4 times */
            while (modificationIndices[residueIndex] == modificaitons.size()) {
                modificationIndices[residueIndex] = -1;
                residueIndex++;
                carryTheZero = true;

                /* exit if we have reached the very end */
                if (residueIndex == modificationIndices.length) {
                    return -1;
                }
            }

            /* if we went from one digit to the next, we increment the most sig.  e.g. 1399 to 1400 */
            if (carryTheZero) {
                modificationIndices[residueIndex]++;
            }

            /* we have reached the end */
            if (modificationIndices[residueIndex] == modificaitons.size()) {
                return -1;
            }

            /* see if we should keep going */
            if (modificationIndices[residueIndex] == -1) {
                keepGoing = true;
            } else {
                /* keep going if this mod has nothing to do with this amino acid */
                keepGoing = (sequence[residueIndex] != modificaitons.get(modificationIndices[residueIndex]).getAminoAcid());
            }


        }

        /* if we went from one digit to the next, we've got to go back to the least significant.  999 goes to 1000, but then that's back to 1001,1002 and so on */
        if (carryTheZero) {
            residueIndex = 0;
        }

        return residueIndex;

    }


}
