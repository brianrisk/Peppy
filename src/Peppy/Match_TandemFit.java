package Peppy;

import Math.MathFunctions;

public class Match_TandemFit extends Match{
	
	
	public String getScoringMethodName() {return "TandemFit";}

	public void calculateScore() {
		byte [] acidSequence = peptide.getAcidSequence();

		int i;
		boolean atLeastOneMatch = false;
		double theoreticalPeakMassLeft, peakMass;
		int peakIndex, seqIndex;
		ionMatchTally = 0;
		
		//we want -1 because most of these spectra will have a match with 
		//the last theoretical peak
		int peptideLengthMinusOne = acidSequence.length - 1;
		
		double [] bIonMatchesWithHighestIntensity = new double[peptideLengthMinusOne];
		double [] yIonMatchesWithHighestIntensity = new double[peptideLengthMinusOne];

		//find the ranges around our theoretical peptides where we
		//count spectrum peaks
		double [] theoreticalPeaksLeft = new double[peptideLengthMinusOne];
		double [] theoreticalPeaksRight = new double[peptideLengthMinusOne];
		
		
		/* y-ion  */
		//computing the left and right boundaries for the ranges where our peaks should land
		theoreticalPeakMassLeft = peptide.getMass() + Properties.rightIonDifference - Properties.peakDifferenceThreshold;
		double peakDifferenceThresholdArea = Properties.peakDifferenceThreshold + Properties.peakDifferenceThreshold;
		for (i = 0; i < peptideLengthMinusOne; i++) {
			theoreticalPeakMassLeft -= AminoAcids.getWeightMono(acidSequence[i]);
			theoreticalPeaksLeft[i] = theoreticalPeakMassLeft;
			theoreticalPeaksRight[i] = theoreticalPeakMassLeft + peakDifferenceThresholdArea;
		}
		
		peakIndex = spectrum.getPeakCount() - 1;
		seqIndex = 0;
		while (peakIndex >= 0) {
			peakMass = spectrum.getPeak(peakIndex).getMass();
			spectrum.getPeak(peakIndex).used = false;
			while (peakMass < theoreticalPeaksLeft[seqIndex]) {
				seqIndex++;
				if (seqIndex == peptideLengthMinusOne) break;
			}
			if (seqIndex == peptideLengthMinusOne) break;
			if (peakMass >= theoreticalPeaksLeft[seqIndex] && peakMass <= theoreticalPeaksRight[seqIndex]) {
				spectrum.getPeak(peakIndex).used = true;
				atLeastOneMatch = true;
				if (yIonMatchesWithHighestIntensity[seqIndex] < spectrum.getPeak(peakIndex).getIntensity()) {
					yIonMatchesWithHighestIntensity[seqIndex] = spectrum.getPeak(peakIndex).getIntensity();
				}
			}
			
			peakIndex--;
		}

		//if 0 matches so far, just get out.
		if (!atLeastOneMatch) {
			score = 0.0;
			return;
		}
			
		
		/* b-ion  */
		theoreticalPeakMassLeft = Properties.leftIonDifference  - Properties.peakDifferenceThreshold;
		for (i = 0; i < peptideLengthMinusOne; i++) {
			theoreticalPeakMassLeft += AminoAcids.getWeightMono(acidSequence[i]);
			theoreticalPeaksLeft[i] = theoreticalPeakMassLeft;
			theoreticalPeaksRight[i] = theoreticalPeakMassLeft + peakDifferenceThresholdArea;
		}
		
		peakIndex = 0;
		seqIndex = 0;
		while (peakIndex < spectrum.getPeakCount()) {
			if (!spectrum.getPeak(peakIndex).used) {
				peakMass = spectrum.getPeak(peakIndex).getMass();
				while (peakMass > theoreticalPeaksRight[seqIndex]) {
					seqIndex++;
					if (seqIndex == peptideLengthMinusOne) break;
				}
				if (seqIndex == peptideLengthMinusOne) break;
				if (peakMass >= theoreticalPeaksLeft[seqIndex] && peakMass <= theoreticalPeaksRight[seqIndex]) {
					atLeastOneMatch = true;
					if (bIonMatchesWithHighestIntensity[seqIndex] < spectrum.getPeak(peakIndex).getIntensity()) {
						bIonMatchesWithHighestIntensity[seqIndex] = spectrum.getPeak(peakIndex).getIntensity();
					}
				}
			}
			peakIndex++;
		}
		
		//if 0 matches so far, just get out.
		if (!atLeastOneMatch) {
			score = 0.0;
			return;
		}
		
		//find out final tally
		boolean yIonTrue, bIonTrue;
		double siblingSum = 0, noSiblingSum = 0;
		for (i = 0; i < peptideLengthMinusOne; i++) {
			
			yIonTrue = yIonMatchesWithHighestIntensity[i] > 0.0;
			bIonTrue = bIonMatchesWithHighestIntensity[i] > 0.0;
			if (yIonTrue) {
				if (bIonTrue) {
					siblingSum += Math.pow(yIonMatchesWithHighestIntensity[i], Properties.peakIntensityExponent);
				} else {
					noSiblingSum += Math.pow(yIonMatchesWithHighestIntensity[i], Properties.peakIntensityExponent);
				}
				ionMatchTally++;
			}
			if (bIonTrue) {
				if (yIonTrue) {
					siblingSum += Math.pow(bIonMatchesWithHighestIntensity[i], Properties.peakIntensityExponent);
				} else {
					noSiblingSum += Math.pow(bIonMatchesWithHighestIntensity[i], Properties.peakIntensityExponent);
				}
				ionMatchTally++;
			}
		}
		score = 1.1 * siblingSum + 0.9 * noSiblingSum;
		
		//This slightly punishes for peptides with many amino acids as they have a slightly higher
		//probability of getting alignments by chance
		score /= MathFunctions.cachedLog(acidSequence.length);
	}

}
