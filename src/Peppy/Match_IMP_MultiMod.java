package Peppy;

import java.util.ArrayList;

import Math.MathFunctions;

public class Match_IMP_MultiMod extends Match {
	
	double bestIMP;
	Modification [] bestModifications;
	ModificationPossibilities modificationPossibilities;

	public Match_IMP_MultiMod(Spectrum spectrum, Peptide peptide, ModificationPossibilities modificationPossibilities) {
		super(spectrum, peptide, false);
		this.modificationPossibilities = modificationPossibilities;
	}
	
	public void calculateScore() {
		findBestIMP(Definitions.WATER_MONO, 0, 0, new Modification[peptide.getLength()]);
		impValue = bestIMP;
		score = -Math.log(bestIMP);
	}
	
	public String getScoringMethodName() {return "IMP MultiMod";}
	
	/**
	 * A recursive functions which, as it is exploring, compares the IMP values found at the
	 * leaf nodes to the bestIMP and, if it is better, updates the bestIMP as well as the
	 * modification array which produced that IMP
	 * @param massTotal
	 * @param modificationMassTotal
	 * @param index
	 * @param modifications
	 */
	private void findBestIMP(double massTotal, double modificationMassTotal, int index, Modification [] modifications) {
		// get the list of possible modifications for the amino acid at this index
		ArrayList<Modification> acidModifications = getModificationsForAminoAcid(peptide.getAcidSequence()[index]);
		
		//variables
		int newIndex = index + 1;
		double imp;
		double newMassTotal;
		double newModificationMassTotal;
		double newMassWithoutMod = massTotal + AminoAcids.getWeightMono(peptide.getAcidSequence()[index]);
		
		
		for (Modification modification: acidModifications) {
			modifications[index] = modification;
			newMassTotal = newMassWithoutMod + modification.getMonoMass();
			newModificationMassTotal = modificationMassTotal + modification.getMonoMass();
			if (newMassTotal > spectrum.getMassPlusMargin()) {
				return;
			} else {
				/*
				 * if we have reached the end of our recursive chain, then
				 * calculate all IMPs and set the best IMP and the modification array
				 * if necessary.
				 */
				if (index == peptide.getLengthMinusOne()) {
					//since this is the final amino acid, the mass sum should be close to our precursor
					if (newMassTotal < spectrum.getMassMinusMargin()) {
						return;
					}
					
					imp = calculateIMP(modifications, newModificationMassTotal);
					//if IMP is less than 0, that means it is an invalid match
					if (imp > 0) {
						if (imp < bestIMP) {
							bestIMP = imp;
							//save the modifications
							for (int i = 0; i < modifications.length; i++) {
								bestModifications[i] = modifications[i];
							}
						}
					}
				}
				/*
				 * If we are still in the midst of our recursion
				 */
				else {
					findBestIMP(newMassTotal, newModificationMassTotal, newIndex, modifications);
				}
			}
		}
	}
	
	private ArrayList<Modification> getModificationsForAminoAcid(byte aminoAcid) {
		return modificationPossibilities.getModificationList(aminoAcid);
	}
	
	
	public double calculateIMP(Modification [] modifications, double modificationMassTotal) {
		byte [] acidSequence = peptide.getAcidSequence();

		int i;
		boolean atLeastOneMatch = false;
		double theoreticalPeakMass, peakMass;
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
		theoreticalPeakMass = modificationMassTotal + peptide.getMass() + Properties.rightIonDifference;
		for (i = 0; i < peptideLengthMinusOne; i++) {
			theoreticalPeakMass -= AminoAcids.getWeightMono(acidSequence[i]);
			theoreticalPeakMass -= modifications[i].getMonoMass();
			theoreticalPeaksLeft[i] = theoreticalPeakMass - Properties.peakDifferenceThreshold;
			theoreticalPeaksRight[i] = theoreticalPeakMass + Properties.peakDifferenceThreshold;	
		}
		
		peakIndex = spectrum.getPeakCount() - 1;
		seqIndex = 0;
		while (peakIndex >= 0) {
			peakMass = spectrum.getPeak(peakIndex).getMass();
			while (peakMass < theoreticalPeaksLeft[seqIndex]) {
				seqIndex++;
				if (seqIndex == peptideLengthMinusOne) break;
			}
			if (seqIndex == peptideLengthMinusOne) break;
			if (peakMass >= theoreticalPeaksLeft[seqIndex] && peakMass <= theoreticalPeaksRight[seqIndex]) {
				atLeastOneMatch = true;
				if (yIonMatchesWithHighestIntensity[seqIndex] < spectrum.getPeak(peakIndex).getIntensity()) {
					yIonMatchesWithHighestIntensity[seqIndex] = spectrum.getPeak(peakIndex).getIntensity();
				}
			}
			
			peakIndex--;
		}

		//if 0 matches so far, just get out.
		if (!atLeastOneMatch) {
			impValue = 1;
			return impValue;
		}
			
		
		/* b-ion  */
		theoreticalPeakMass = Properties.leftIonDifference;
		for (i = 0; i < peptideLengthMinusOne; i++) {
			theoreticalPeakMass += AminoAcids.getWeightMono(acidSequence[i]);
			theoreticalPeakMass += modifications[i].getMonoMass();
			theoreticalPeaksLeft[i] = theoreticalPeakMass - Properties.peakDifferenceThreshold;
			theoreticalPeaksRight[i] = theoreticalPeakMass + Properties.peakDifferenceThreshold;
		}
		
		peakIndex = 0;
		seqIndex = 0;
		while (peakIndex < spectrum.getPeakCount()) {
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
			
			peakIndex++;
		}
		
		//if 0 matches so far, just get out.
		if (!atLeastOneMatch) {
			impValue = 1;
			return impValue;
		}
		
		//find out ionMatchTally and intensity distributions
		int totalIonsAbove50 = 0;
		int totalIonsAbove25 = 0;
		int totalIonsAbove12 = 0;
		int totalIonsAbove06 = 0;
		double totalMatchingIntensity = 0.0;
		boolean yIonMatch, bIonMatch, ionMatch, previousIonMatch = false;
		int acidMatchTally = 0;
		int consecutiveAcidTally = 0;
		for (i = 0; i < peptideLengthMinusOne; i++) {
			yIonMatch = yIonMatchesWithHighestIntensity[i] > 0.0;
			bIonMatch = bIonMatchesWithHighestIntensity[i] > 0.0;
			ionMatch = bIonMatch || yIonMatch;
			if (ionMatch) acidMatchTally++;
			if (previousIonMatch && ionMatch) consecutiveAcidTally++;
			if (yIonMatch) {
				ionMatchTally++;
				totalMatchingIntensity += yIonMatchesWithHighestIntensity[i];
				if (yIonMatchesWithHighestIntensity[i] > spectrum.getMedianIntensity()) {
					totalIonsAbove50++;
					if (yIonMatchesWithHighestIntensity[i] > spectrum.getIntensity25Percent()) {
						totalIonsAbove25++;
						if (yIonMatchesWithHighestIntensity[i] > spectrum.getIntensity12Percent()) {
							totalIonsAbove12++;
							if (yIonMatchesWithHighestIntensity[i] > spectrum.getIntensity06Percent()) {
								totalIonsAbove06++;
							}
						}
					}
				}
			}
			if (bIonMatch) {
				ionMatchTally++;
				totalMatchingIntensity += bIonMatchesWithHighestIntensity[i];
				if (bIonMatchesWithHighestIntensity[i] > spectrum.getMedianIntensity()) {
					totalIonsAbove50++;
					if (bIonMatchesWithHighestIntensity[i] > spectrum.getIntensity25Percent()) {
						totalIonsAbove25++;
						if (bIonMatchesWithHighestIntensity[i] > spectrum.getIntensity12Percent()) {
							totalIonsAbove12++;
							if (bIonMatchesWithHighestIntensity[i] > spectrum.getIntensity06Percent()) {
								totalIonsAbove06++;
							}
						}
					}
				}
			}
			previousIonMatch = ionMatch;
			
		}
		
		
		//Variables for binomial probabilities
		int n;
		int k;
		double p;
		
		//peak match probability is the binomial distribution
		n = acidSequence.length * 2;
		k = ionMatchTally;
		p = spectrum.getCoverage();
		double peakMatchProbability = MathFunctions.getBinomialProbability(n, k, p);
		
		
		if (peakMatchProbability > 0.5) {
			impValue = 1;
			return impValue;
		}
		
		//consecutive match probability
		n = peptideLengthMinusOne;
		k = consecutiveAcidTally;
		p = (double) acidMatchTally / peptideLengthMinusOne;
		double consecutiveProbability = MathFunctions.getBinomialProbability(n, k, p);
		
		
		//probability of ions being above thresholds
		double intensityProbability = 1;
		if (ionMatchTally > 0) {
			n = ionMatchTally;
			k = totalIonsAbove50;
			intensityProbability = MathFunctions.getCachedBinomialProbability50(n, k);
			if (totalIonsAbove50 > 0) {
				n = totalIonsAbove50;
				k = totalIonsAbove25;
				intensityProbability *= MathFunctions.getCachedBinomialProbability50(n, k);
				if (totalIonsAbove25 > 0) {
					n = totalIonsAbove25;
					k = totalIonsAbove12;
					intensityProbability *= MathFunctions.getCachedBinomialProbability50(n, k);
					if (totalIonsAbove12 > 0) {
						n = totalIonsAbove12;
						k = totalIonsAbove06;
						intensityProbability *= MathFunctions.getCachedBinomialProbability50(n, k);
					}
				}
			}
		}
		
		impValue =  peakMatchProbability  * intensityProbability * consecutiveProbability;
		return impValue;
	}

}
