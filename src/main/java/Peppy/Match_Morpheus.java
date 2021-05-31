package Peppy;

/**
 * Implementation of the Morpheus PSM score
 * 
 * 
 * 
 * Copyright 2013, Brian Risk
 * 
 * @author Brian Risk
 *
 */
public class Match_Morpheus extends Match {

	@Override
	public void calculateScore() {
		score = 0;
		byte [] acidSequence = peptide.getAcidSequence();
		
		int peptideLengthMinusOne = acidSequence.length - 1;
		if (acidSequence[peptideLengthMinusOne] == AminoAcids.STOP) peptideLengthMinusOne--;
		
		//will hold our peak boundaries
		double [] theoreticalPeaksLeft = new double[peptideLengthMinusOne];
		double [] theoreticalPeaksRight = new double[peptideLengthMinusOne];
		
		//find y ions
		double [] yIonMatchesWithHighestIntensity = findYIons(peptideLengthMinusOne, theoreticalPeaksLeft, theoreticalPeaksRight);
		
		//find b ions
		double [] bIonMatchesWithHighestIntensity = findBIons(peptideLengthMinusOne, theoreticalPeaksLeft, theoreticalPeaksRight);
		
		double totalMatchedIntensity = 0;
		for (double intensity: yIonMatchesWithHighestIntensity) {
			if (intensity > 0) {
				score += 1;
				totalMatchedIntensity += intensity;
			}
		}
		for (double intensity: bIonMatchesWithHighestIntensity) {
			if (intensity > 0) {
				score += 1;
				totalMatchedIntensity += intensity;
			}
		}
		
		/* adding the matched intensity percentage */
		score += totalMatchedIntensity / matchesSpectrum.getSpectrum().getTotalIntensity();
		
	}

}
