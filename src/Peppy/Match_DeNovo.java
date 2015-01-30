package Peppy;

import Math.MassError;

/**
 * Experimental class for DeNovo sequencing
 * 
 * @author brianrisk
 *
 */
public class Match_DeNovo extends Match {
	
	Peptide peptide;
	Spectrum spectrum;
	double nTerminalMass;
	double cTerminalMass;

	public Match_DeNovo(Peptide peptide, Spectrum spectrum, double nTerminalMass, double cTerminalMass) {
		this.peptide = peptide;
		this.spectrum = spectrum;
		this.nTerminalMass = nTerminalMass;
		this.cTerminalMass = cTerminalMass;
	}
	
	@Override
	public void calculateScore() {

		int i;
		boolean atLeastOneMatch = false;
		double theoreticalPeakMass, peakMass;
		int peakIndex, seqIndex;
		
		//we want -1 because most of these spectra will have a match with 
		//the last theoretical peak
		int peptideLengthMinusOne = peptide.getLengthMinusOne();
		
		double [] bIonMatchesWithHighestIntensity = new double[peptideLengthMinusOne];
		double [] yIonMatchesWithHighestIntensity = new double[peptideLengthMinusOne];

		//find the ranges around our theoretical peptides where we
		//count spectrum peaks
		double [] theoreticalPeaksLeft = new double[peptideLengthMinusOne];
		double [] theoreticalPeaksRight = new double[peptideLengthMinusOne];
		
		
		/* y-ion  */
		//computing the left and right boundaries for the ranges where our peaks should land
		theoreticalPeakMass = nTerminalMass + cTerminalMass + peptide.getMass() + Properties.rightIonDifference;
		if (Properties.isITRAQ) {
			theoreticalPeakMass -= Properties.ITRAQ_REAGENT;
		}
		for (i = 0; i < peptideLengthMinusOne; i++) {
			theoreticalPeakMass -= peptide.getResidueMass(i);
			if (i == 0) theoreticalPeakMass -= cTerminalMass;
			if (i == peptideLengthMinusOne) theoreticalPeakMass -= nTerminalMass;
			theoreticalPeaksLeft[i] = theoreticalPeakMass - MassError.getDaltonError(Properties.fragmentTolerance, theoreticalPeakMass);
			theoreticalPeaksRight[i] = theoreticalPeakMass + MassError.getDaltonError(Properties.fragmentTolerance, theoreticalPeakMass);	
		}
		
		peakIndex = matchesSpectrum.getSpectrum().getPeakCount() - 1;
		seqIndex = 0;
		while (peakIndex >= 0) {
			peakMass = matchesSpectrum.getSpectrum().getPeak(peakIndex).getMass();
			while (peakMass < theoreticalPeaksLeft[seqIndex]) {
				seqIndex++;
				if (seqIndex == peptideLengthMinusOne) break;
			}
			if (seqIndex == peptideLengthMinusOne) break;
			if (peakMass >= theoreticalPeaksLeft[seqIndex] && peakMass <= theoreticalPeaksRight[seqIndex]) {
				atLeastOneMatch = true;
				if (yIonMatchesWithHighestIntensity[seqIndex] < matchesSpectrum.getSpectrum().getPeak(peakIndex).getIntensity()) {
					yIonMatchesWithHighestIntensity[seqIndex] = matchesSpectrum.getSpectrum().getPeak(peakIndex).getIntensity();
				}
			}
			
			peakIndex--;
		}

		//if 0 matches so far, just get out.
		if (!atLeastOneMatch) {
			score = 0;
		}
			
		
		/* b-ion  */
		theoreticalPeakMass = Properties.leftIonDifference;
		if (Properties.isITRAQ) {
			theoreticalPeakMass += Properties.ITRAQ_REAGENT;
		}
		for (i = 0; i < peptideLengthMinusOne; i++) {
			theoreticalPeakMass += peptide.getResidueMass(i);
			if (i == 0) theoreticalPeakMass += cTerminalMass;
			if (i == peptideLengthMinusOne) theoreticalPeakMass += nTerminalMass;
			theoreticalPeaksLeft[i] = theoreticalPeakMass - MassError.getDaltonError(Properties.fragmentTolerance, theoreticalPeakMass);
			theoreticalPeaksRight[i] = theoreticalPeakMass + MassError.getDaltonError(Properties.fragmentTolerance, theoreticalPeakMass);
		}
		
		peakIndex = 0;
		seqIndex = 0;
		while (peakIndex < matchesSpectrum.getSpectrum().getPeakCount()) {
			peakMass = matchesSpectrum.getSpectrum().getPeak(peakIndex).getMass();
			while (peakMass > theoreticalPeaksRight[seqIndex]) {
				seqIndex++;
				if (seqIndex == peptideLengthMinusOne) break;
			}
			if (seqIndex == peptideLengthMinusOne) break;
			if (peakMass >= theoreticalPeaksLeft[seqIndex] && peakMass <= theoreticalPeaksRight[seqIndex]) {
				atLeastOneMatch = true;
				if (bIonMatchesWithHighestIntensity[seqIndex] < matchesSpectrum.getSpectrum().getPeak(peakIndex).getIntensity()) {
					bIonMatchesWithHighestIntensity[seqIndex] = matchesSpectrum.getSpectrum().getPeak(peakIndex).getIntensity();
				}
			}
			
			peakIndex++;
		}
		
		
//		return calculateIMP(peptideLengthMinusOne, yIonMatchesWithHighestIntensity, bIonMatchesWithHighestIntensity);
	}

}
