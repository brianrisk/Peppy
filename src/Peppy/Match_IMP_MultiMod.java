package Peppy;

import java.util.ArrayList;

public class Match_IMP_MultiMod extends Match {
	
	double bestIMP = 1;
	Modification [] bestModifications;
	
	public void calculateScore() {
		bestModifications = new Modification[peptide.getLength()];
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
		ArrayList<Modification> acidModifications = Properties.modificationPossibilities.getModificationList(peptide.getAcidSequence()[index]);
		
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
			
			//if theoretical peptide mass is too heavy, exit
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
					
					if (imp < bestIMP) {
						bestIMP = imp;
						//save the modifications
						for (int i = 0; i < modifications.length; i++) {
							bestModifications[i] = modifications[i];
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
	

	/**
	 * This eventually calls Match's calculateIMP, but first calculates b and y ions
	 * given the list of modifications.
	 * 
	 * @param modifications
	 * @param modificationMassTotal
	 * @return
	 */
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
			theoreticalPeaksLeft[i] = theoreticalPeakMass - Properties.fragmentTolerance;
			theoreticalPeaksRight[i] = theoreticalPeakMass + Properties.fragmentTolerance;	
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
			theoreticalPeaksLeft[i] = theoreticalPeakMass - Properties.fragmentTolerance;
			theoreticalPeaksRight[i] = theoreticalPeakMass + Properties.fragmentTolerance;
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
		
		return calculateIMP(peptideLengthMinusOne, yIonMatchesWithHighestIntensity, bIonMatchesWithHighestIntensity);
	}

	public Modification [] getModifications() {
		return bestModifications;
	}
	

}
