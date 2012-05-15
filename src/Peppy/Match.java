package Peppy;

import Math.MassError;
import Math.MathFunctions;

/**
 * An abstract object which can be extended to implement scoring mechanisms
 * 
 * IMP is natively implemented with this class as it is also used as a validity
 * check for E values
 * @author Brian Risk
 *
 */
public abstract class Match implements Comparable<Match>{
	
	protected MatchesSpectrum matchesSpectrum;
	protected Peptide peptide;
	
	protected double score = 0.0;

	protected double impValue = -1;
	
	/* I'm setting up "cScore" to be a way to help sort out the situation where 
	 * matching 25 out of 50 AAs in a peptide produces a better score than matching
	 * 10 out of 10 AAs in a shorter peptide 
	 * 
	 * This could, potentially, be a method to eliminate the false positives produced
	 * by modified matches outscoring the modified that have equal quality of match
	 * 
	 * */
	protected double cScore = -1;
	
	protected int numberOfIonsMatched = 0;
	
	/* so that we can combine sets of matches but still be able to separate them out */
	private int trackingIdentifier = 0;
	
	private static int sortTracker = 0;
	public final static int SORT_BY_SCORE = sortTracker++;
	public final static int SORT_BY_SPECTRUM_ID_THEN_SCORE = sortTracker++;
	public final static int SORT_BY_LOCUS = sortTracker++;
	public final static int SORT_BY_SPECTRUM_PEPTIDE_MASS_DIFFERENCE = sortTracker++;
	public final static int SORT_BY_SPECTRUM_ID_THEN_PEPTIDE = sortTracker++;
	public final static int SORT_BY_IMP_VALUE = sortTracker++;
	public final static int SORT_BY_PEPTIDE_THEN_SCORE = sortTracker++;
	
	//default is that we sort matches by score
	private static int sortParameter = SORT_BY_SCORE;
	
	

	
	public abstract void calculateScore();
	public abstract String getScoringMethodName();
	
	public void recalculateIMP() {
		impValue = -1;
		calculateIMP();
	}
	
	public double calculateIMP() {
		if (impValue < 0) {
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
				impValue = 1;
				return impValue;
			}
			
			//find b ions
			double [] bIonMatchesWithHighestIntensity = findBIons(peptideLengthMinusOne, theoreticalPeaksLeft, theoreticalPeaksRight);
			
			calculateIMP(peptideLengthMinusOne, yIonMatchesWithHighestIntensity, bIonMatchesWithHighestIntensity);
			
		}
		return impValue;
	}
	
	public double getIMP() {
		if (impValue < 0) {
			calculateIMP();
		}
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
		int ionMatchTally = 0;
		for (int i = 0; i < peptideLengthMinusOne; i++) {
			yIonMatch = yIonMatchesWithHighestIntensity[i] > 0.0;
			bIonMatch = bIonMatchesWithHighestIntensity[i] > 0.0;
			if (yIonMatchesWithHighestIntensity[i] > bIonMatchesWithHighestIntensity[i]) appropriateIonIsMoreIntenseTally++;
			if (yIonMatch) {
				ionMatchTally++;
//				totalMatchingIntensity += yIonMatchesWithHighestIntensity[i];
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
//				totalMatchingIntensity += bIonMatchesWithHighestIntensity[i];
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
		
		//Variables for binomial probabilities
		int n;
		int k;
		double p;
		
		/* peak match probability is the binomial distribution */
		p = matchesSpectrum.getSpectrum().getCoverage();
		double targetMass = matchesSpectrum.getSpectrum().getMass() - AminoAcids.getWeightMono(peptide.getAcidSequence()[0]) - AminoAcids.getWeightMono(peptide.getAcidSequence()[peptideLengthMinusOne]);
//		p /= targetMass;
		p /= matchesSpectrum.getSpectrum().getMass();
		
		/* sometimes p will be unrealistically small due to high accuracy fragment masses */
//		double minimumP = Definitions.PROTON_MONO / targetMass;
//		if (p < minimumP) p = minimumP;
		
		n = peptideLengthMinusOne * 2;
		k = ionMatchTally;
		double peakMatchProbability = MathFunctions.getBinomialProbability(n, k, p);
//		double peakMatchProbability = Math.pow(p, k);
		
		/* setting our "coverage" score */
		if (k > numberOfIonsMatched) numberOfIonsMatched = k;
		
		//TODO: optimize this.  it is a great place to improve performance
		//NOTE:  if we end up changing the binomial probability to return the log, this will need to get changed
		if (peakMatchProbability > 0.25) {
			impValue = 1;
			return 1;
		}
		
		//y greater than b probability
		n = peptideLengthMinusOne;
		k = appropriateIonIsMoreIntenseTally;
		double appropriateIonIsMoreIntenseProbablity = MathFunctions.getCachedBinomialProbability50(n, k);
		
		
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

		

		//finding the probability of the total of the matching peak intensities
//		int totalPeakCombnations = 100;
//		double rise = MathFunctions.cachedLog(totalPeakCombnations);
//		double maxScore = spectrum.getMaxValueForCombination(ionMatchTally);
//		double minScore = ionMatchTally * spectrum.getMinimumIntensity();
//		double run = (maxScore - minScore);
//		double totalMatchingIntensityProbability = -1 * rise  * (totalMatchingIntensity - minScore) / run + rise;
//		totalMatchingIntensityProbability = Math.exp(totalMatchingIntensityProbability - totalPeakCombnations);
//		U.p(totalMatchingIntensityProbability);
		
		impValue =  peakMatchProbability  * intensityProbability * appropriateIonIsMoreIntenseProbablity;
		
		
		if (impValue > 1) impValue = 1;
		return impValue;
	}
	
	
	
	protected double [] findYIons(int peptideLengthMinusOne, double [] theoreticalPeaksLeft, double [] theoreticalPeaksRight) {
		byte [] acidSequence = peptide.getAcidSequence();
		
		int i;
		boolean atLeastOneMatch = false;
		double theoreticalPeakMass, peakMass;
		int peakIndex, seqIndex;
		
		double [] yIonMatchesWithHighestIntensity = new double[peptideLengthMinusOne];
		
		/* y-ion  */
		//computing the left and right boundaries for the ranges where our peaks should land
		theoreticalPeakMass = peptide.getMass() + Properties.rightIonDifference;
		for (i = 0; i < peptideLengthMinusOne; i++) {
			theoreticalPeakMass -= AminoAcids.getWeightMono(acidSequence[i]);
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

	protected double [] findBIons(int peptideLengthMinusOne, double [] theoreticalPeaksLeft, double [] theoreticalPeaksRight) {
		byte [] acidSequence = peptide.getAcidSequence();
		
		int i;
		double theoreticalPeakMass, peakMass;
		int peakIndex, seqIndex;
		
		double [] bIonMatchesWithHighestIntensity = new double[peptideLengthMinusOne];
		
		/* b-ion  */
		theoreticalPeakMass = Properties.leftIonDifference;
		for (i = 0; i < peptideLengthMinusOne; i++) {
			theoreticalPeakMass += AminoAcids.getWeightMono(acidSequence[i]);
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
			}
			peakIndex++;
		}
		
		return bIonMatchesWithHighestIntensity;
	}

	public int compareTo(Match match) {
			if (sortParameter == SORT_BY_SCORE) {
				//want to sort from greatest to least
				if (score > match.getScore()) return -1;
				if (score < match.getScore()) return  1;
				return 0;
			} else
			if (sortParameter == SORT_BY_LOCUS) {
				if (peptide.getParentSequence().getId() < match.getPeptide().getParentSequence().getId()) return -1;
				if (peptide.getParentSequence().getId() > match.getPeptide().getParentSequence().getId()) return  1;
				//in case sequences equal, compare index
				if(peptide.getStartIndex() < match.getPeptide().getStartIndex()) return -1;
				if(peptide.getStartIndex() > match.getPeptide().getStartIndex()) return  1;
				return 0;
			} else
			if (sortParameter == SORT_BY_IMP_VALUE) {
				if (impValue < match.getIMP()) return -1;
				if (impValue > match.getIMP()) return  1;
				return 0;
			} else
			if (sortParameter == SORT_BY_SPECTRUM_ID_THEN_SCORE) {
				if (matchesSpectrum.getSpectrum().getId() < match.getSpectrum().getId()) return -1;
				if (matchesSpectrum.getSpectrum().getId() > match.getSpectrum().getId()) return  1;
				//if spectrum is sorted, also sort by tandemFit
				if (score > match.getScore()) return -1;
				if (score < match.getScore()) return  1;
				return 0;
			} else
			if (sortParameter == SORT_BY_SPECTRUM_ID_THEN_PEPTIDE) {
				//first by spectrum ID
				if (matchesSpectrum.getSpectrum().getId() < match.getSpectrum().getId()) return -1;
				if (matchesSpectrum.getSpectrum().getId() > match.getSpectrum().getId()) return  1;
				//then by start location
				if(peptide.getStartIndex() < match.getPeptide().getStartIndex()) return -1;
				if(peptide.getStartIndex() > match.getPeptide().getStartIndex()) return  1;
				//then by alphabetical order of peptides
				int shortLength = peptide.getAcidSequence().length;
				if (match.getPeptide().getAcidSequence().length < shortLength) shortLength = match.getPeptide().getAcidSequence().length;
				for (int i = 0; i < shortLength; i++) {
					if (match.getPeptide().getAcidSequence()[i] != peptide.getAcidSequence()[i]) return match.getPeptide().getAcidSequence()[i] - peptide.getAcidSequence()[i];
				}
				return 0;
			} else				
			if (sortParameter == SORT_BY_SPECTRUM_PEPTIDE_MASS_DIFFERENCE) {
				//i'm putting this calculation in here as this is a not-often-used sort
				//so calculating and storing this for every match is unnecessary
				double myDifference = matchesSpectrum.getSpectrum().getMass() - peptide.getMass();
				double theirDifference = match.getSpectrum().getMass() - match.getPeptide().getMass();
				if (myDifference > theirDifference) return -1;
				if (myDifference < theirDifference) return  1;
				return 0;
			} else 
			if (sortParameter == SORT_BY_PEPTIDE_THEN_SCORE) {
				/* first sorting by mass */
				if (getPeptide().getMass() > match.getPeptide().getMass()) return  1;
				if (getPeptide().getMass() < match.getPeptide().getMass()) return -1;
				
				/* sorting alphabetically */
				int shortLength = peptide.getAcidSequence().length;
				if (match.getPeptide().getAcidSequence().length < shortLength) shortLength = match.getPeptide().getAcidSequence().length;
				for (int i = 0; i < shortLength; i++) {
					if (match.getPeptide().getAcidSequence()[i] != peptide.getAcidSequence()[i]) return match.getPeptide().getAcidSequence()[i] - peptide.getAcidSequence()[i];
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
				if (score < match.getScore()) return  1;
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
	 * Returns the default score.  This could be TandemFit or HMM.
	 */
	public double getScore() {
		return score;
	}

	
	public double getCScore() {
		return numberOfIonsMatched * numberOfIonsMatched / peptide.getLengthMinusOne();
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
		sb.append(getSpectrum().getMD5());
		sb.append('\t');
		sb.append(getSpectrum().getFile().getAbsolutePath());
		sb.append('\t');
		sb.append(getScore());
		sb.append('\t');
		sb.append(getSpectrum().getPrecursorMZ());
		sb.append('\t');
		sb.append(getSpectrum().getMass());
		sb.append('\t');
		sb.append(getPeptide().getAcidSequenceString());
		sb.append('\t');
		sb.append(getPeptide().getStartIndex());
		sb.append('\t');
		sb.append(getPeptide().getStopIndex());
		sb.append('\t');
		if (Properties.isSequenceFileDNA) {
			sb.append(U.getFileNameWithoutSuffix(getPeptide().getParentSequence().getSequenceFile()));
			sb.append('\t');
			sb.append(getPeptide().getIntronStartIndex());
			sb.append('\t');
			sb.append(getPeptide().getIntronStopIndex());
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
		sb.append(matchesSpectrum.getMatches().size());
		sb.append('\t');
		sb.append(getSpectrum().getCharge());
		sb.append('\t');
		sb.append(getPeptide().getCleavageAcidCount());
		sb.append('\t');
		sb.append(getPeptide().isInORF());
		sb.append('\t');
		sb.append(getPeptide().getORFSize());
		sb.append('\t');
		sb.append(getPeptide().getHydrophobicProportion());
		sb.append('\t');
		sb.append(getPeptide().getHydrophilicProportion());
		sb.append('\t');
		sb.append(hasModification());
		sb.append('\t');
		sb.append(getMoificationdMass());
		sb.append('\t');
		sb.append(getModificationIndex());
		sb.append('\t');
		sb.append(isModificationLocationCertain());
		sb.append('\t');
		sb.append(getCScore());
		return sb.toString();
	}
	

	public boolean hasModification() {
		return false;
	}
	
	/**
	 * A match might exceed the precursor tolerance of a non-modified search, but still
	 * wouldn't be what we would consider to be modified.  This addresses that distinction.
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
	
	public ModificationEntry [] getModifications() {
		return new ModificationEntry [peptide.getLength()];
	}
	
	public int getNumberOfModifications() {
		int tally = 0;
		for (ModificationEntry mod: getModifications()) {
			if (mod == null) continue;
			if (mod.getMonoMass() == 0) continue;
			tally++;
		}
		return tally;
	}
	
	/**
	 * 
	 * @return peptide mass - spectrum mass
	 */
	public double getMassDifference() {
		return peptide.getMass() - matchesSpectrum.getSpectrum().getMass();
	}
	
	public int getTrackingIdentifier() {
		return trackingIdentifier;
	}

	public void setTrackingIdentifier(int trackingIdentifier) {
		this.trackingIdentifier = trackingIdentifier;
	}
	

}
