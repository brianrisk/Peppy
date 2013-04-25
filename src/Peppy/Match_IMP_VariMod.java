package Peppy;

import Math.MassError;

/**
 * Copyright 2012, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class Match_IMP_VariMod extends Match {
	
	//Peptide/Spectrum difference in mass
	double modificationMass;
	int modificationIndex = 0;
	
	/* this boolean tracks if the location of the modification is ambiguous */
	boolean modificationLocationCertain = true;
	
	
	/* true if mass difference between peptide and spectrum is significant */
	boolean hasModification;
	
	//this is which amino acid we think has the modification
	boolean [] potentialPlacesForModification;
	
	public Match_IMP_VariMod() {}
	

	

	@Override
	public void calculateScore() {
		
		
		/* modification mass is the difference between our precursor and our theoretical mass */
		modificationMass = matchesSpectrum.getSpectrum().getMass() - peptide.getMass();
		
		/* calculate the score */
		score = -Math.log10(calculateIMP());
		/*
		 * if a modification is less than the fragment tolerance, then it would not affect the score when
		 * compared to a non-modified search.
		 */
		hasModification = Math.abs(modificationMass) > MassError.getDaltonError(Properties.fragmentTolerance, matchesSpectrum.getSpectrum().getMass());
		
		
		/* this punishment will serve as a tie breaker for comparing
		 * matches with modifications and those without.  If we find a spectrum matching
		 * two different peptides with the same score, we should believe the unmodified
		 * match more.
		 * 
		 * the if is to avoid less than zero values as this would create a bug
		 * when trying to find our histogram for e values.
		 * we don't really worry about ties at this level of score anyway.
		 */
		if (hasModification) {
			if (score > 3) {
				score -= 3;
			}
		}
		
	}
	
	
	public double calculateIMP() {
		if (impValue < 0) {
			double impValue1 = calculateIMP(modificationMass, 0);
			double impValue2 = calculateIMP(modificationMass, peptide.getAcidSequence().length - 1);
			
			/* if we see some promise, explore further */
			//TODO: come up with a method for justifying 1.0E-8 
			if (impValue1 < 1.0E-8 || impValue2 < 1.0E-8) {
				double bestIMP = impValue1;
				double tempIMP;
				for (int i = 1; i < peptide.getAcidSequence().length - 1; i++) {
					tempIMP = calculateIMP(modificationMass, i);
					if (tempIMP < bestIMP) {
						bestIMP = tempIMP;
						modificationIndex = i;
						modificationLocationCertain = true;
					} else {
						if (tempIMP == bestIMP) {
							modificationLocationCertain = false;
						}
					}
				}
				if (impValue2 < bestIMP) {
					bestIMP = impValue2;
					modificationIndex =  peptide.getAcidSequence().length - 1;
					modificationLocationCertain = true;
				}
				impValue = bestIMP;
				
			/* the score if not initially a good IMP */
			} else {
				if (impValue2 < impValue1) {
					impValue = impValue2;
				} else {
					impValue = impValue1;
				}
			}
		}
		return impValue;
	}
	
	
	/**
	 * This eventually calls Match's calculateIMP, but first calculates b and y ions
	 * given the modification.
	 * 
	 * @param modificationMass
	 * @param modifiedIndex
	 * @return
	 */
	protected double calculateIMP(double modificationMass, int modifiedIndex) {

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
		theoreticalPeakMass = modificationMass + peptide.getMass() + Properties.rightIonDifference;
		if (Properties.isITRAQ) {
			theoreticalPeakMass -= Definitions.ITRAQ_REAGENT;
		}
		for (i = 0; i < peptideLengthMinusOne; i++) {
			theoreticalPeakMass -= peptide.getResidueMass(i);
			if (i == modifiedIndex) theoreticalPeakMass -= modificationMass;
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
			impValue = 1;
			return impValue;
		}
			
		
		/* b-ion  */
		theoreticalPeakMass = Properties.leftIonDifference;
		if (Properties.isITRAQ) {
			theoreticalPeakMass += Definitions.ITRAQ_REAGENT;
		}
		for (i = 0; i < peptideLengthMinusOne; i++) {
			theoreticalPeakMass += peptide.getResidueMass(i);
			if (i == modifiedIndex) theoreticalPeakMass += modificationMass;
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
		
		
		return calculateIMP(peptideLengthMinusOne, yIonMatchesWithHighestIntensity, bIonMatchesWithHighestIntensity);
	}
	
	public double getMoificationdMass() {return modificationMass;}
	
	public int getModificationIndex() {
		return modificationIndex;
	}

	public boolean hasModification() {
		return hasModification;
	}
	
	public boolean isFromModificationSearches() {
		return (Math.abs(modificationMass) > MassError.getDaltonError(Properties.precursorTolerance, matchesSpectrum.getSpectrum().getMass()));
	}
	
	public boolean isModificationLocationCertain() {
		return modificationLocationCertain;
	}




	public int getNumberOfModifications() {
		if (hasModification()) {
			return 1;
		} else {
			return 0;
		}
	}
	
//	public String toString() {
//		String out = super.toString();
//		StringBuffer sb = new StringBuffer(out);
//		if (Properties.searchModifications) {
//			sb.append('\t');
//			sb.append(hasModification());
//			sb.append('\t');
//			sb.append(getMoificationdMass());
//			sb.append('\t');
//			sb.append(getModificationIndex() + 1);
//		}
//		return sb.toString();
//	}

}
