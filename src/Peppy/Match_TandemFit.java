package Peppy;

import Math.MathFunctions;

public class Match_TandemFit extends Match{
	
	
	public String getScoringMethodName() {return "TandemFit";}

	public void calculateScore() {
		/*
		 * PART 1:  Find matching peaks
		 */
		ionMatchTally = 0;
		byte [] acidSequence = peptide.getAcidSequence();
		
		//we want -1 because most of these spectra will have a match with 
		//the last theoretical peak
		int peptideLengthMinusOne = acidSequence.length - 1;
		if (acidSequence[peptideLengthMinusOne] == AminoAcids.STOP) peptideLengthMinusOne--;
		
		//will hold our peak boundaries
		double [] theoreticalPeaksLeft = new double[peptideLengthMinusOne];
		double [] theoreticalPeaksRight = new double[peptideLengthMinusOne];
		
		//find y ions
		double [] yIonMatchesWithHighestIntensity = findYIons(peptideLengthMinusOne, theoreticalPeaksLeft, theoreticalPeaksRight);
		
		//if null is returned, that means no y ions match, which is automatic disqualification
		if (yIonMatchesWithHighestIntensity == null) {
			score = 0;
			return;
		}
		
		//find b ions
		double [] bIonMatchesWithHighestIntensity = findBIons(peptideLengthMinusOne, theoreticalPeaksLeft, theoreticalPeaksRight);
		
		/*
		 * PART 2: Apply TandemFit formula
		 */
		boolean yIonTrue, bIonTrue;
		double siblingSum = 0, noSiblingSum = 0;
		for (int i = 0; i < peptideLengthMinusOne; i++) {
			
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
