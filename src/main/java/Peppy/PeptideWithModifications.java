package Peppy;


/**
 * Copyright 2013, Brian Risk
 *
 * @author Brian Risk
 */
public class PeptideWithModifications extends Peptide {

    private double[] modifications;

    public PeptideWithModifications(String sequence, double[] modifications) {
        super(sequence);
        this.modifications = modifications;
        for (double mod : modifications) {
            mass += mod;
        }
    }

    public double[] getModifications() {
        return modifications;
    }

    public double getResidueMass(int index) {
        return AminoAcids.getWeightMono(acidSequence[index]) + modifications[index];
    }

}
