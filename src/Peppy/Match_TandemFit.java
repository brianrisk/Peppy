package Peppy;

import Math.MathFunctions;

public class Match_TandemFit extends Match{
	
	public Match_TandemFit(Spectrum spectrum, Peptide peptide) {
		super(spectrum, peptide);
	}
	
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
		}
		
		//find out final tally
		boolean yIonTrue, bIonTrue;
		double amountToAdd;
		for (i = 0; i < peptideLengthMinusOne; i++) {
			
			yIonTrue = yIonMatchesWithHighestIntensity[i] > 0.0;
			bIonTrue = bIonMatchesWithHighestIntensity[i] > 0.0;
			if (yIonTrue) {
//				score += Math.pow(yIonMatchesWithHighestIntensity[i], Properties.peakIntensityExponent);
				if (bIonTrue) {
//					amountToAdd =  1.17 * Math.pow(yIonMatchesWithHighestIntensity[i], Properties.peakIntensityExponent);
					amountToAdd =  Properties.YBtrue * Math.pow(yIonMatchesWithHighestIntensity[i], Properties.peakIntensityExponent);
				} else {
//					amountToAdd =  1.43 * Math.pow(yIonMatchesWithHighestIntensity[i], Properties.peakIntensityExponent);
					amountToAdd =  Properties.YBfalse * Math.pow(yIonMatchesWithHighestIntensity[i], Properties.peakIntensityExponent);
				}
//				amountToAdd *=  1.3;
				score += amountToAdd;
				ionMatchTally++;
			}
			if (bIonTrue) {
				//count b-ions less if not y-ion present, more if present
				if (yIonTrue) {
//					amountToAdd =  1.1 * Math.pow(bIonMatchesWithHighestIntensity[i], Properties.peakIntensityExponent);
					amountToAdd =  Properties.BYtrue * Math.pow(bIonMatchesWithHighestIntensity[i], Properties.peakIntensityExponent);
				} else {
//					amountToAdd =  0.9 * Math.pow(bIonMatchesWithHighestIntensity[i], Properties.peakIntensityExponent);
					amountToAdd =  Properties.BYfalse * Math.pow(bIonMatchesWithHighestIntensity[i], Properties.peakIntensityExponent);
				}
				score += amountToAdd;
				ionMatchTally++;
			}
		}
		//This slightly punishes for peptides with many amino acids as they have a slightly higher
		//probability of getting alignments by chance
		score /= MathFunctions.cachedLog(acidSequence.length);
	}

}
