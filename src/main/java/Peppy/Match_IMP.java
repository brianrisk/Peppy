package Peppy;

import Experimental.PaperExperiments;
import Math.MathFunctions;


/**
 * 
 * "IMP" stands for "ion match probability"
 * 
 * The algorithm is describe in:
 * "A peptide-spectrum scoring system based on ion alignment, intensity and pair probabilities." 
 * Risk BA, Edwards NJ, Giddings MC. Journal of proteome research (July 22, 2013).
 * 
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class Match_IMP extends Match {


	@Override
	public void calculateScore() {
		calculateIMP();	
	}
	
	public double calculateIMP() {
		double impValue;
		
		byte [] acidSequence = peptide.getAcidSequence();
		
		int peptideLengthMinusOne = acidSequence.length - 1;
		if (acidSequence[peptideLengthMinusOne] == AminoAcids.STOP) peptideLengthMinusOne--;
		
		//will hold our peak boundaries
		double [] theoreticalPeaksLeft = new double[peptideLengthMinusOne];
		double [] theoreticalPeaksRight = new double[peptideLengthMinusOne];
		
		//find y ions
		double [] yIonMatchesWithHighestIntensity = findYIons(peptideLengthMinusOne, theoreticalPeaksLeft, theoreticalPeaksRight);
		
		//if null is returned, that means no y ions match, which is automatic disqualification
		if (yIonMatchesWithHighestIntensity == null) {
			impValue = 1;
			return impValue;
		}
		
		//find b ions
		double [] bIonMatchesWithHighestIntensity = findBIons(peptideLengthMinusOne, theoreticalPeaksLeft, theoreticalPeaksRight);
		
		impValue = calculateIMP(peptideLengthMinusOne, yIonMatchesWithHighestIntensity, bIonMatchesWithHighestIntensity);
			
		return impValue;
	}
	

	
	public void setSpectrumMatches(MatchesSpectrum matchesSpectrum) {
		this.matchesSpectrum = matchesSpectrum;
	}

	public void setPeptide(Peptide peptide) {
		this.peptide = peptide;
	}

	/**
	 * Here our y and b ion matches have been found so what remains is to calculate the IMP
	 * 
	 * @param peptideLengthMinusOne
	 * @param yIonMatchesWithHighestIntensity
	 * @param bIonMatchesWithHighestIntensity
	 */
	protected double calculateIMP(int peptideLengthMinusOne, double [] yIonMatchesWithHighestIntensity, double [] bIonMatchesWithHighestIntensity) {
		//find out ionMatchTally and intensity distributions
		int totalIonsAbove50 = 0;
		int totalIonsAbove25 = 0;
		int totalIonsAbove12 = 0;
		int totalIonsAbove06 = 0;
		
		/* the totalMatchingIntensity variably MIGHT be very 
		 * useful in the future as an extra probability to incorporate 
		 * don't get rid of it just yet! */
//		double totalMatchingIntensity = 0.0;
		boolean yIonMatch, bIonMatch = false;
		int appropriateIonIsMoreIntenseTally = 0;
		int pairCount = 0;
		int ionMatchTally = 0;
		for (int i = 0; i < peptideLengthMinusOne; i++) {
			yIonMatch = yIonMatchesWithHighestIntensity[i] > 0.0;
			bIonMatch = bIonMatchesWithHighestIntensity[i] > 0.0;
			
			if (yIonMatch || bIonMatch) pairCount++;
			if (yIonMatch) {
				
				if (yIonMatchesWithHighestIntensity[i] >= bIonMatchesWithHighestIntensity[i]) appropriateIonIsMoreIntenseTally++;
				ionMatchTally++;
				if (yIonMatchesWithHighestIntensity[i] > matchesSpectrum.getSpectrum().getMedianIntensity()) {
					totalIonsAbove50++;
					if (yIonMatchesWithHighestIntensity[i] > matchesSpectrum.getSpectrum().getIntensity25Percent()) {
						totalIonsAbove25++;
						if (yIonMatchesWithHighestIntensity[i] > matchesSpectrum.getSpectrum().getIntensity12Percent()) {
							totalIonsAbove12++;
							if (yIonMatchesWithHighestIntensity[i] > matchesSpectrum.getSpectrum().getIntensity06Percent()) {
								totalIonsAbove06++;
							}
						}
					}
				}
			}
			if (bIonMatch) {
				ionMatchTally++;
				if (bIonMatchesWithHighestIntensity[i] > matchesSpectrum.getSpectrum().getMedianIntensity()) {
					totalIonsAbove50++;
					if (bIonMatchesWithHighestIntensity[i] > matchesSpectrum.getSpectrum().getIntensity25Percent()) {
						totalIonsAbove25++;
						if (bIonMatchesWithHighestIntensity[i] > matchesSpectrum.getSpectrum().getIntensity12Percent()) {
							totalIonsAbove12++;
							if (bIonMatchesWithHighestIntensity[i] > matchesSpectrum.getSpectrum().getIntensity06Percent()) {
								totalIonsAbove06++;
							}
						}
					}
				}
			}
		}
		/*
		 * If no ions matched, that's a bad match
		 */
		if (ionMatchTally <= 0) {
			return 1;
		}
		
		// Variables for binomial probabilities
		int n;
		int k;
		double p;
		
		/* peak match probability is the binomial distribution */
		p = matchesSpectrum.getSpectrum().getCoverage();		
		n = peptideLengthMinusOne * 2;
		k = ionMatchTally;
		double peakMatchProbability = MathFunctions.getBinomialProbability(n, k, p);
		
		// NOTE:  if we end up changing the binomial probability to return the log, this will need to get changed
		if (peakMatchProbability > 0.01) {
			return 1;
		}
		
		// y greater than b probability
		n = pairCount;
		k = appropriateIonIsMoreIntenseTally;
		double appropriateIonIsMoreIntenseProbablity = MathFunctions.getCachedBinomialProbability50(n, k);

		
		
		// probability of ions being above thresholds
		double intensityProbability = 1;
		if (totalIonsAbove50 > 0) {
			n = ionMatchTally;
//			n = peptideLengthMinusOne * 2;
			k = totalIonsAbove50;
			if (PaperExperiments.p2a)
			intensityProbability = MathFunctions.getCachedBinomialProbability50(n, k);
			if (totalIonsAbove25 > 0) {
				n = totalIonsAbove50;
				k = totalIonsAbove25;
				if (PaperExperiments.p2b)
				intensityProbability *= MathFunctions.getCachedBinomialProbability50(n, k);
				if (totalIonsAbove12 > 0) {
					n = totalIonsAbove25;
					k = totalIonsAbove12;
					if (PaperExperiments.p2c)
					intensityProbability *= MathFunctions.getCachedBinomialProbability50(n, k);
					if (totalIonsAbove06 > 0) {
						n = totalIonsAbove12;
						k = totalIonsAbove06;
						if (PaperExperiments.p2d)
						intensityProbability *= MathFunctions.getCachedBinomialProbability50(n, k);
					}
				}
			}
		}


		double impValue  = peakMatchProbability * intensityProbability * appropriateIonIsMoreIntenseProbablity;

		if (impValue > 1) impValue = 1;
		
		score = -Math.log10(impValue);
		
		return impValue;
	}
	
	


}
