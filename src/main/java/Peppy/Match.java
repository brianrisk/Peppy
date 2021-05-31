package Peppy;

import Math.MassError;

import java.util.ArrayList;

/**
 * An abstract object which can be extended to implement scoring mechanisms
 * <p>
 * IMP is natively implemented with this class as it is also used as a validity
 * check for E values
 * <p>
 * Copyright 2013, Brian Risk
 *
 * @author Brian Risk
 */
public abstract class Match implements Comparable<Match> {

    protected MatchesSpectrum matchesSpectrum;
    protected Peptide peptide;

    protected double score = 0.0;

    private static int sortTracker = 0;
    public final static int SORT_BY_SCORE = sortTracker++;
    public final static int SORT_BY_SPECTRUM_ID_THEN_SCORE = sortTracker++;
    public final static int SORT_BY_LOCUS = sortTracker++;
    public final static int SORT_BY_SPECTRUM_PEPTIDE_MASS_DIFFERENCE = sortTracker++;
    public final static int SORT_BY_SPECTRUM_ID_THEN_PEPTIDE = sortTracker++;
    public final static int SORT_BY_PEPTIDE_THEN_SCORE = sortTracker++;

    //default is that we sort matches by score
    private static int sortParameter = SORT_BY_SCORE;


    public abstract void calculateScore();


    protected double[] findYIons(int peptideLengthMinusOne, double[] theoreticalPeaksLeft, double[] theoreticalPeaksRight) {

        int i;
        boolean atLeastOneMatch = false;
        double theoreticalPeakMass, peakMass;
        int peakIndex, seqIndex;

        double[] yIonMatchesWithHighestIntensity = new double[peptideLengthMinusOne];

        /* y-ion  */
        //computing the left and right boundaries for the ranges where our peaks should land
        theoreticalPeakMass = peptide.getMass() + Properties.rightIonDifference + Properties.modNTerminal + Properties.modCTerminal;
        if (Properties.isITRAQ) {
            theoreticalPeakMass -= Properties.ITRAQ_REAGENT;
        }
        for (i = 0; i < peptideLengthMinusOne; i++) {
            theoreticalPeakMass -= peptide.getResidueMass(i);
            theoreticalPeaksLeft[i] = theoreticalPeakMass - MassError.getDaltonError(Properties.fragmentTolerance, theoreticalPeakMass);
            theoreticalPeaksRight[i] = theoreticalPeakMass + MassError.getDaltonError(Properties.fragmentTolerance, theoreticalPeakMass);
        }

        peakIndex = matchesSpectrum.getSpectrum().getPeakCount() - 1;
        seqIndex = 0;
        while (peakIndex >= 0) {
            peakMass = matchesSpectrum.getSpectrum().getPeak(peakIndex).getMass();
            matchesSpectrum.getSpectrum().getPeak(peakIndex).used = false;
            while (peakMass < theoreticalPeaksLeft[seqIndex]) {
                seqIndex++;
                if (seqIndex == peptideLengthMinusOne) break;
            }
            if (seqIndex == peptideLengthMinusOne) break;
            if (peakMass >= theoreticalPeaksLeft[seqIndex] && peakMass <= theoreticalPeaksRight[seqIndex]) {
                matchesSpectrum.getSpectrum().getPeak(peakIndex).used = true;
                atLeastOneMatch = true;
                if (yIonMatchesWithHighestIntensity[seqIndex] < matchesSpectrum.getSpectrum().getPeak(peakIndex).getIntensity()) {
                    yIonMatchesWithHighestIntensity[seqIndex] = matchesSpectrum.getSpectrum().getPeak(peakIndex).getIntensity();
                }
            }

            peakIndex--;
        }

        //if 0 matches so far, just get out.
        if (!atLeastOneMatch) {
            return null;
        }

        return yIonMatchesWithHighestIntensity;
    }


    protected double[] findBIons(int peptideLengthMinusOne, double[] theoreticalPeaksLeft, double[] theoreticalPeaksRight) {

        int i;
        double theoreticalPeakMass, peakMass;
        int peakIndex, seqIndex;

        double[] bIonMatchesWithHighestIntensity = new double[peptideLengthMinusOne];

        /* b-ion  */
        theoreticalPeakMass = Properties.leftIonDifference + Properties.modNTerminal;
        if (Properties.isITRAQ) {
            theoreticalPeakMass += Properties.ITRAQ_REAGENT;
        }
        for (i = 0; i < peptideLengthMinusOne; i++) {
            theoreticalPeakMass += peptide.getResidueMass(i);
            theoreticalPeaksLeft[i] = theoreticalPeakMass - MassError.getDaltonError(Properties.fragmentTolerance, theoreticalPeakMass);
            theoreticalPeaksRight[i] = theoreticalPeakMass + MassError.getDaltonError(Properties.fragmentTolerance, theoreticalPeakMass);
        }

        peakIndex = 0;
        seqIndex = 0;
        while (peakIndex < matchesSpectrum.getSpectrum().getPeakCount()) {
            if (!matchesSpectrum.getSpectrum().getPeak(peakIndex).used) {
                peakMass = matchesSpectrum.getSpectrum().getPeak(peakIndex).getMass();
                while (peakMass > theoreticalPeaksRight[seqIndex]) {
                    seqIndex++;
                    if (seqIndex == peptideLengthMinusOne) break;
                }
                if (seqIndex == peptideLengthMinusOne) break;
                if (peakMass >= theoreticalPeaksLeft[seqIndex] && peakMass <= theoreticalPeaksRight[seqIndex]) {
                    if (bIonMatchesWithHighestIntensity[seqIndex] < matchesSpectrum.getSpectrum().getPeak(peakIndex).getIntensity()) {
                        bIonMatchesWithHighestIntensity[seqIndex] = matchesSpectrum.getSpectrum().getPeak(peakIndex).getIntensity();
                    }
                }
            } else {
                /* as we marked peaks finding y-ions, we must now unmark them */
                matchesSpectrum.getSpectrum().getPeak(peakIndex).used = false;
            }
            peakIndex++;
        }


        return bIonMatchesWithHighestIntensity;
    }

    public ArrayList<Double> getFragmentErrors() {

        ArrayList<Double> out = new ArrayList<Double>();

        byte[] acidSequence = peptide.getAcidSequence();

        int peptideLengthMinusOne = acidSequence.length - 1;
        if (acidSequence[peptideLengthMinusOne] == AminoAcids.STOP) peptideLengthMinusOne--;

        //will hold our peak boundaries
        double[] theoreticalPeaksLeft = new double[peptideLengthMinusOne];
        double[] theoreticalPeaksRight = new double[peptideLengthMinusOne];

        int i;
        double theoreticalPeakMass, peakMass;
        int peakIndex, seqIndex;

        double[] yIonMatchesWithHighestIntensity = new double[peptideLengthMinusOne];
        double[] theoreticalPeakMasses = new double[peptideLengthMinusOne];

        /* y-ion  */
        //computing the left and right boundaries for the ranges where our peaks should land
        theoreticalPeakMass = peptide.getMass() + Properties.rightIonDifference + Properties.modNTerminal + Properties.modCTerminal;
        if (Properties.isITRAQ) {
            theoreticalPeakMass -= Properties.ITRAQ_REAGENT;
        }
        for (i = 0; i < peptideLengthMinusOne; i++) {
            theoreticalPeakMass -= peptide.getResidueMass(i);
            theoreticalPeakMasses[i] = theoreticalPeakMass;
            theoreticalPeaksLeft[i] = theoreticalPeakMass - MassError.getDaltonError(Properties.fragmentTolerance, theoreticalPeakMass);
            theoreticalPeaksRight[i] = theoreticalPeakMass + MassError.getDaltonError(Properties.fragmentTolerance, theoreticalPeakMass);
        }

        peakIndex = matchesSpectrum.getSpectrum().getPeakCount() - 1;
        seqIndex = 0;
        double fragmentError = 0;
        while (peakIndex >= 0) {
            peakMass = matchesSpectrum.getSpectrum().getPeak(peakIndex).getMass();
            while (peakMass < theoreticalPeaksLeft[seqIndex]) {
                seqIndex++;
                if (fragmentError != 0) {
                    out.add(fragmentError);
                    fragmentError = 0;
                }
                if (seqIndex == peptideLengthMinusOne) break;
            }
            if (seqIndex == peptideLengthMinusOne) break;
            if (peakMass >= theoreticalPeaksLeft[seqIndex] && peakMass <= theoreticalPeaksRight[seqIndex]) {
                if (yIonMatchesWithHighestIntensity[seqIndex] < matchesSpectrum.getSpectrum().getPeak(peakIndex).getIntensity()) {
                    yIonMatchesWithHighestIntensity[seqIndex] = matchesSpectrum.getSpectrum().getPeak(peakIndex).getIntensity();

                    /* calculate fragment error */
                    fragmentError = MassError.getPPMDifference(theoreticalPeakMasses[seqIndex], peakMass);
                }
            }

            peakIndex--;
        }

        if (fragmentError != 0) {
            out.add(fragmentError);
        }

        return out;


    }


    protected int getNumberOfIonsMatched() {
        byte[] acidSequence = peptide.getAcidSequence();

        int peptideLengthMinusOne = acidSequence.length - 1;
        if (acidSequence[peptideLengthMinusOne] == AminoAcids.STOP) peptideLengthMinusOne--;

        //will hold our peak boundaries
        double[] theoreticalPeaksLeft = new double[peptideLengthMinusOne];
        double[] theoreticalPeaksRight = new double[peptideLengthMinusOne];

        //find y ions
        double[] yIonMatchesWithHighestIntensity = findYIons(peptideLengthMinusOne, theoreticalPeaksLeft, theoreticalPeaksRight);


        //find b ions
        double[] bIonMatchesWithHighestIntensity = findBIons(peptideLengthMinusOne, theoreticalPeaksLeft, theoreticalPeaksRight);

        return getNumberOfIonsMatched(peptideLengthMinusOne, yIonMatchesWithHighestIntensity, bIonMatchesWithHighestIntensity);
    }

    protected int getNumberOfIonsMatched(int peptideLengthMinusOne, double[] yIonMatchesWithHighestIntensity, double[] bIonMatchesWithHighestIntensity) {
        int ionMatchTally = 0;
        for (int i = 0; i < peptideLengthMinusOne; i++) {
            if (yIonMatchesWithHighestIntensity[i] > 0.0) ionMatchTally++;
            if (bIonMatchesWithHighestIntensity[i] > 0.0) ionMatchTally++;
        }
        return ionMatchTally;
    }


    public int compareTo(Match match) {
        if (sortParameter == SORT_BY_SCORE) {
            //want to sort from greatest to least
            if (score > match.getScore()) return -1;
            if (score < match.getScore()) return 1;
            return 0;
        } else if (sortParameter == SORT_BY_LOCUS) {
            if (peptide.getParentSequence().getId() < match.getPeptide().getParentSequence().getId()) return -1;
            if (peptide.getParentSequence().getId() > match.getPeptide().getParentSequence().getId()) return 1;
            //in case sequences equal, compare index
            if (peptide.getStartIndex() < match.getPeptide().getStartIndex()) return -1;
            if (peptide.getStartIndex() > match.getPeptide().getStartIndex()) return 1;
            return 0;
        } else if (sortParameter == SORT_BY_SPECTRUM_ID_THEN_SCORE) {
            if (matchesSpectrum.getSpectrum().getId() < match.getSpectrum().getId()) return -1;
            if (matchesSpectrum.getSpectrum().getId() > match.getSpectrum().getId()) return 1;
            //if spectrum is sorted, also sort by tandemFit
            if (score > match.getScore()) return -1;
            if (score < match.getScore()) return 1;
            return 0;
        } else if (sortParameter == SORT_BY_SPECTRUM_ID_THEN_PEPTIDE) {
            //first by spectrum ID
            if (matchesSpectrum.getSpectrum().getId() < match.getSpectrum().getId()) return -1;
            if (matchesSpectrum.getSpectrum().getId() > match.getSpectrum().getId()) return 1;
            //then by start location
            if (peptide.getStartIndex() < match.getPeptide().getStartIndex()) return -1;
            if (peptide.getStartIndex() > match.getPeptide().getStartIndex()) return 1;
            //then by alphabetical order of peptides
            int shortLength = peptide.getAcidSequence().length;
            if (match.getPeptide().getAcidSequence().length < shortLength)
                shortLength = match.getPeptide().getAcidSequence().length;
            for (int i = 0; i < shortLength; i++) {
                if (match.getPeptide().getAcidSequence()[i] != peptide.getAcidSequence()[i])
                    return match.getPeptide().getAcidSequence()[i] - peptide.getAcidSequence()[i];
            }
            return 0;
        } else if (sortParameter == SORT_BY_SPECTRUM_PEPTIDE_MASS_DIFFERENCE) {
            //i'm putting this calculation in here as this is a not-often-used sort
            //so calculating and storing this for every match is unnecessary
            double myDifference = matchesSpectrum.getSpectrum().getMass() - peptide.getMass();
            double theirDifference = match.getSpectrum().getMass() - match.getPeptide().getMass();
            if (myDifference > theirDifference) return -1;
            if (myDifference < theirDifference) return 1;
            return 0;
        } else if (sortParameter == SORT_BY_PEPTIDE_THEN_SCORE) {
            /* first sorting by mass */
            if (getPeptide().getMass() > match.getPeptide().getMass()) return 1;
            if (getPeptide().getMass() < match.getPeptide().getMass()) return -1;

            /* sorting alphabetically */
            int shortLength = peptide.getAcidSequence().length;
            if (match.getPeptide().getAcidSequence().length < shortLength)
                shortLength = match.getPeptide().getAcidSequence().length;
            for (int i = 0; i < shortLength; i++) {
                if (match.getPeptide().getAcidSequence()[i] != peptide.getAcidSequence()[i])
                    return match.getPeptide().getAcidSequence()[i] - peptide.getAcidSequence()[i];
            }

            /* sorting by score */
            if (score < match.getScore()) return 1;
            if (score > match.getScore()) return -1;
            return 0;
        } else

        {
            //we want to sort from greatest to least great
            //so -1 is returned where 1 usually is
            if (score > match.getScore()) return -1;
            if (score < match.getScore()) return 1;
            return 0;
        }
    }

    public boolean equals(Match match) {
//		if (getScore() == match.getScore()) {
        if (peptide.getStartIndex() == match.getPeptide().getStartIndex())
            if (peptide.equals(match.getPeptide()))
                if (matchesSpectrum.getSpectrum().equals(match.getSpectrumMatches().getSpectrum()))
                    return true;
//		}
        return false;
    }

    /**
     * Returns the score.
     */
    public double getScore() {
        return score;
    }


    /**
     * @return the spectrum
     */
    public Spectrum getSpectrum() {
        return matchesSpectrum.getSpectrum();
    }

    public MatchesSpectrum getSpectrumMatches() {
        return matchesSpectrum;
    }


    /**
     * @return the peptide
     */
    public Peptide getPeptide() {
        return peptide;
    }


    public void setScore(double score) {
        this.score = score;
    }


    public static void setSortParameter(int sortParameter) {
        Match.sortParameter = sortParameter;
    }

    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append(getSpectrum().getId());
        sb.append('\t');
        sb.append(getSpectrum().getFileLocus());
        sb.append('\t');
//		sb.append(getSpectrum().getMD5());
//		sb.append('\t');
        sb.append(getSpectrum().getFile().getAbsolutePath());
        sb.append('\t');
        sb.append(getScore());
        sb.append('\t');
        sb.append(getPeptide().getMass());
        sb.append('\t');
        sb.append(getSpectrum().getMass());
        sb.append('\t');
        sb.append(getPeptide().getAcidSequenceString());
        sb.append('\t');
        sb.append(peptide.getPreviousAminoAcid());
        sb.append('\t');
        sb.append(getPeptide().getStartIndex());
        sb.append('\t');
        sb.append(getPeptide().getStopIndex());
        sb.append('\t');
        if (Properties.isSequenceFileNucleotide) {
            if (getPeptide().getParentSequence() != null) {
                if (getPeptide().getProtein() != null) {
                    sb.append(getPeptide().getProtein().getName());
                } else {
                    sb.append(U.getFileNameWithoutSuffix(getPeptide().getParentSequence().getSequenceFile()));
                }
            } else {
                if (getPeptide().getProtein() != null) {
                    sb.append(getPeptide().getProtein().getName());
                } else {
                    sb.append("null");
                }
            }
            sb.append('\t');
            sb.append(getPeptide().isForward() ? "+" : "-");
            sb.append('\t');
            sb.append(getPeptide().isSpliced());
        } else {
            if (getPeptide().getProtein() != null) {
                sb.append(getPeptide().getProtein().getName());
            }
        }
        sb.append('\t');
        sb.append(getSpectrum().getCharge());
        sb.append('\t');
        sb.append(hasModification());
        sb.append('\t');
        sb.append(getMoificationdMass());
        sb.append('\t');
        sb.append(getModificationIndex());
        sb.append('\t');
        sb.append(isModificationLocationCertain());
        return sb.toString();
    }

    public boolean hasModification() {
        return false;
    }

    /**
     * A match might exceed the precursor tolerance of a non-modified search, but still
     * wouldn't be what we would consider to be modified.  This addresses that distinction.
     *
     * @return
     */
    public boolean isFromModificationSearches() {
        return false;
    }

    public double getMoificationdMass() {
        return 0;
    }

    public int getModificationIndex() {
        return 0;
    }

    public boolean isModificationLocationCertain() {
        return true;
    }


    /**
     * @return peptide mass - spectrum mass
     */
    public double getMassDifference() {
        return peptide.getMass() - matchesSpectrum.getSpectrum().getMass();
    }


    public void setMatchesSpectrum(MatchesSpectrum matchesSpectrum) {
        this.matchesSpectrum = matchesSpectrum;
    }


    public void setPeptide(Peptide peptide) {
        this.peptide = peptide;
    }


}
