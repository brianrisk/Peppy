package Peppy;

import Math.HasValue;


/**
 * Is a data class that stores:
 * 1) a sequence of amino acids.
 * 2) the theoretical mass
 * 3) the begin index (end index can be calculated)
 * 4) boolean forward (false = reverse)
 * 5) the reading frame
 * <p>
 * Copyright 2013, Brian Risk
 *
 * @author Brian Risk
 */
public class Peptide implements Comparable<Peptide>, HasValue {

    protected byte[] acidSequence;
    protected double mass;
    private int startIndex;
    private int stopIndex;
    private int intronStartIndex;
    private int intronStopIndex;
    private boolean forward;
    private Sequence parentSequence;
    private Protein protein;
    private boolean isSpliced;
    private int lengthMinusOne;
    private boolean inORF;
    private int ORFSize;
    private char previousAminoAcid;

    /* for tracking FDR */
    private boolean isDecoy = false;

    /* for tracking which result set this peptide belongs to */
    private int trackingIdentifier = 0;


    /**
     * just gets an amino acid sequence.
     *
     * @param sequence
     */
    public Peptide(String sequence) {
        this(sequence, 0, sequence.length(), -1, -1, true, null, null, false, false, -1, '.');
    }

    /**
     * This constructor is here so that we can make peptide
     * objects that are just mass.  This allows us to use the
     * Collections binary search to locate where a peptide of
     * a given mass would be in a list of sorted peptides
     *
     * @param mass
     */
    public Peptide(double mass) {
        this.mass = mass;
    }

    public Peptide(String acidSequence, int startIndex, int stopIndex, int intronStartIndex, int intronStopIndex, boolean forward, Sequence parentSequence, boolean isSpliced) {
        this(acidSequence, startIndex, stopIndex, intronStartIndex, intronStopIndex, forward, parentSequence, null, isSpliced, false, -1, '.');
    }

    public Peptide(String acidSequence, int startIndex, int stopIndex, int intronStartIndex, int intronStopIndex, boolean forward, Sequence parentSequence, Protein protein, boolean isSpliced, boolean inORF, int ORFSize, char previousAminoAcid) {
        this.acidSequence = AminoAcids.getByteArrayForString(acidSequence);
        this.mass = calculateMass();
        this.startIndex = startIndex;
        this.stopIndex = stopIndex;
        this.intronStartIndex = intronStartIndex;
        this.intronStopIndex = intronStopIndex;
        this.forward = forward;
        this.parentSequence = parentSequence;
        this.protein = protein;
        this.isSpliced = isSpliced;
        this.lengthMinusOne = this.acidSequence.length - 1;
        this.inORF = inORF;
        this.ORFSize = ORFSize;
        this.previousAminoAcid = previousAminoAcid;
    }


    public int getLengthMinusOne() {
        return lengthMinusOne;
    }

    @Override
    public String toString() {
//		return mass + "\t" + getAcidSequenceString() + "\t" + startIndex + "\t" + proteinName;
        return getAcidSequenceString() + "\t" + getMass() + "\t" + getStartIndex() + "\t" + getStopIndex() + "\t" + forward;
//		return  getAcidSequenceString();
    }


    public int compareTo(Peptide other) {
        if (mass > other.getMass()) return 1;
        if (mass < other.getMass()) return -1;
        return 0;
    }

    /**
     * Okay, this equals is not in line with the way things work for compareTo.
     * this compares acid sequences for equality.  compareTo compares masses.
     * <p>
     * the real trick for equality is ignoring any trailing stop (".") codon
     */
    public boolean equals(byte[] otherAcidSequence) {
        int ourLength = acidSequence.length;
        int theirLength = otherAcidSequence.length;

        //ignoring terminating STOPs
        if (acidSequence[acidSequence.length - 1] == AminoAcids.STOP) ourLength--;
        if (otherAcidSequence[otherAcidSequence.length - 1] == AminoAcids.STOP) theirLength--;

        //if peptids not same length, they are not equal
        if (ourLength != theirLength) return false;

        //compare each acid
        for (int i = 0; i < ourLength; i++) {
            if (acidSequence[i] != otherAcidSequence[i]) return false;
        }

        //if reached this point, they are equal
        return true;
    }

    /**
     * Equals if every sequential acid weighs the same as that of the other sequence
     */
    public boolean equalsByAcidMasses(byte[] otherAcidSequence) {
        if (!equals(otherAcidSequence)) return false;

        for (int i = 0; i < acidSequence.length; i++) {
            if (AminoAcids.getWeightMono(acidSequence[i]) != AminoAcids.getWeightMono(otherAcidSequence[i]))
                return false;
        }

        return true;

    }

    public boolean equals(Peptide peptide) {
        if (mass == peptide.getMass()) {
            return equals(peptide.getAcidSequence());
        } else {
            return false;
        }
    }

    public boolean equals(String acidSequenceString) {
        return equals(AminoAcids.getByteArrayForString(acidSequenceString));
    }


    /**
     * @return the sequence
     */
    public byte[] getAcidSequence() {
        return acidSequence;
    }

    public String getAcidSequenceString() {
        return AminoAcids.getStringForByteArray(acidSequence);
    }

    public double getResidueMass(int index) {
        return AminoAcids.getWeightMono(acidSequence[index]);
    }

    public int getLength() {
        return acidSequence.length;
    }


    /**
     * @return the mass
     */
    public double getMass() {
        return mass;
    }


    /**
     * @return the index
     */
    public int getStartIndex() {
        if (forward) {
            return startIndex;
        } else {
            return stopIndex + 3;
        }
    }

    public int getStopIndex() {
        if (forward) {
            return stopIndex;
        } else {
            return startIndex + 3;
        }
    }

    public int getIntronStartIndex() {
        if (forward) {
            return intronStartIndex;
        } else {
            return intronStopIndex + 1;
        }
    }

    public int getIntronStopIndex() {
        if (forward) {
            return intronStopIndex;
        } else {
            return intronStartIndex + 1;
        }
    }

    public Protein getProtein() {
        return protein;
    }


    /**
     * @return the forward
     */
    public boolean isForward() {
        return forward;
    }


    public boolean isSpliced() {
        return isSpliced;
    }

    public boolean isInORF() {
        return inORF;
    }

    public int getORFSize() {
        return ORFSize;
    }

    public char getPreviousAminoAcid() {
        return previousAminoAcid;
    }

    /**
     * @return the parentSequence
     */
    public Sequence getParentSequence() {
        return parentSequence;
    }

    /**
     * This will calculate either the mono mass or the average mass depending on the setting
     * in your Properties object.
     *
     * @return
     */
    private double calculateMass() {
        double mass = 0.0;
        if (Properties.useMonoMass) {
            for (int i = 0; i < acidSequence.length; i++) {
                if (AminoAcids.isValid(acidSequence[i])) {
                    mass += AminoAcids.getWeightMono(acidSequence[i]);
                } else {
                    mass = -1;
                    return mass;
                }
            }
            mass += Definitions.WATER_MONO;

            /* add reporter ion masses */
            if (Properties.isITRAQ) {
                mass += Properties.ITRAQ_REAGENT;
            }

        } else {
            for (int i = 0; i < acidSequence.length; i++) {
                if (AminoAcids.isValid(acidSequence[i])) {
                    mass += AminoAcids.getWeightAverage(acidSequence[i]);
                } else {
                    mass = -1;
                    return mass;
                }
            }
            mass += Definitions.WATER_AVERAGE;

            /* add reporter ion masses */
            if (Properties.isITRAQ) {
                mass += Properties.ITRAQ_REAGENT;
            }
        }
        return mass;
    }

    public double getValue() {
        return getMass();
    }

    public double getHydrophobicProportion() {
        double out = 0;
        for (int i = 0; i < acidSequence.length; i++) {
            /* (leucine, valine, isoleucine, phenylalanine, methionine, cysteine and tryptophan */
            if (acidSequence[i] == AminoAcids.G) out++;
            if (acidSequence[i] == AminoAcids.A) out++;
            if (acidSequence[i] == AminoAcids.V) out++;
            if (acidSequence[i] == AminoAcids.L) out++;
            if (acidSequence[i] == AminoAcids.I) out++;
            if (acidSequence[i] == AminoAcids.M) out++;
            if (acidSequence[i] == AminoAcids.F) out++;
            if (acidSequence[i] == AminoAcids.W) out++;
            if (acidSequence[i] == AminoAcids.P) out++;
        }
        out /= acidSequence.length;
        return out;
    }

    public double getHydrophilicProportion() {
        double out = 0;
        for (int i = 0; i < acidSequence.length; i++) {
            if (acidSequence[i] == AminoAcids.S) out++;
            if (acidSequence[i] == AminoAcids.T) out++;
            if (acidSequence[i] == AminoAcids.C) out++;
            if (acidSequence[i] == AminoAcids.Y) out++;
            if (acidSequence[i] == AminoAcids.N) out++;
            if (acidSequence[i] == AminoAcids.Q) out++;
        }
        out /= acidSequence.length;
        return out;
    }


    public boolean isDecoy() {
        return isDecoy;
    }

    public void setDecoy(boolean isDecoy) {
        this.isDecoy = isDecoy;
    }

    public int getTrackingIdentifier() {
        return trackingIdentifier;
    }

    public void setTrackingIdentifier(int trackingIdentifier) {
        this.trackingIdentifier = trackingIdentifier;
    }


}
