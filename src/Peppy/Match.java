package Peppy;

import java.util.ArrayList;

import Math.HasEValue;
import Math.MathFunctions;
import Utilities.U;

/**
 * An abstract object which can be extended to implement scoring mechanisms
 * 
 * IMP is natively implemented with this class as it is also used as a validity
 * check for E values
 * @author Brian Risk
 *
 */
public abstract class Match implements Comparable<Match>, HasEValue{
	
	protected Spectrum spectrum;
	protected Peptide peptide;
	
	protected double score = 0.0;
	protected int ionMatchTally = 0;
	
	public double scoreRatio = -1;
	public int repeatCount = 0; 
	public int rank = Integer.MAX_VALUE;
	
	private double eValue;
	protected double impValue = -1;
	
	private static int sortTracker = 0;
	public final static int SORT_BY_SCORE = sortTracker++;
	public final static int SORT_BY_SPECTRUM_ID_THEN_SCORE = sortTracker++;
	public final static int SORT_BY_LOCUS = sortTracker++;
	public final static int SORT_BY_SCORE_RATIO = sortTracker++;
	public final static int SORT_BY_E_VALUE = sortTracker++;
	public final static int SORT_BY_RANK_THEN_E_VALUE = sortTracker++;
	public final static int SORT_BY_SPECTRUM_PEPTIDE_MASS_DIFFERENCE = sortTracker++;
	public final static int SORT_BY_SPECTRUM_ID_THEN_PEPTIDE = sortTracker++;
	public final static int SORT_BY_RANK_THEN_SCORE = sortTracker++;
	public final static int SORT_BY_IMP_VALUE = sortTracker++;
	
	//default is that we sort matches by score
	private static int sortParameter = SORT_BY_SCORE;
	
	public Match(Spectrum spectrum, Peptide peptide) {
		this(spectrum, peptide, true);
	}
	
	public Match(Spectrum spectrum, Peptide peptide, boolean calculateScore) {
		this.spectrum = spectrum;
		this.peptide = peptide;
		if (calculateScore) calculateScore();
	}
	
	public abstract void calculateScore();
	public abstract String getScoringMethodName();
	
	public static Match createMatch(Spectrum spectrum, Peptide peptide) {
		if (Properties.defaultScore == Properties.DEFAULT_SCORE_IMP) {
			return new Match_IMP(spectrum, peptide);
		} else {
			if (Properties.defaultScore == Properties.DEFAULT_SCORE_TANDEM_FIT) {
				return new Match_TandemFit(spectrum, peptide);
			} else {
				if (Properties.defaultScore == Properties.DEFAULT_SCORE_HMM) {
					return new Match_HMM(spectrum, peptide);
				} else {
					U.p("ERROR: undefined scoring method");
					System.exit(1);
				}
			}
		}
		return null;
	}

	
	
	public double calculateIMP() {
		if (impValue < 0) {
			ionMatchTally = 0;
			
//			double [] bIonMatchesWithHighestIntensity = findBIons();
			
			
			byte [] acidSequence = peptide.getAcidSequence();
	
			int i;
			boolean atLeastOneMatch = false;
			double theoreticalPeakMass, peakMass;
			int peakIndex, seqIndex;
			
			
			//we want -1 because most of these spectra will have a match with 
			//the last theoretical peak
			int peptideLengthMinusOne = acidSequence.length - 1;
			if (acidSequence[peptideLengthMinusOne] == AminoAcids.STOP) peptideLengthMinusOne--;
			
			double [] bIonMatchesWithHighestIntensity = new double[peptideLengthMinusOne];
			double [] yIonMatchesWithHighestIntensity = new double[peptideLengthMinusOne];
	
			//find the ranges around our theoretical peptides where we
			//count spectrum peaks
			double [] theoreticalPeaksLeft = new double[peptideLengthMinusOne];
			double [] theoreticalPeaksRight = new double[peptideLengthMinusOne];
			
			
			/* y-ion  */
			//computing the left and right boundaries for the ranges where our peaks should land
			theoreticalPeakMass = peptide.getMass() + Properties.rightIonDifference;
			for (i = 0; i < peptideLengthMinusOne; i++) {
				theoreticalPeakMass -= AminoAcids.getWeightMono(acidSequence[i]);
				theoreticalPeaksLeft[i] = theoreticalPeakMass - Properties.peakDifferenceThreshold;
				theoreticalPeaksRight[i] = theoreticalPeakMass + Properties.peakDifferenceThreshold;
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
				impValue = 1;
				return impValue;
			}
				
			
			/* b-ion  */
			theoreticalPeakMass = Properties.leftIonDifference;
			for (i = 0; i < peptideLengthMinusOne; i++) {
				theoreticalPeakMass += AminoAcids.getWeightMono(acidSequence[i]);
				theoreticalPeaksLeft[i] = theoreticalPeakMass - Properties.peakDifferenceThreshold;
				theoreticalPeaksRight[i] = theoreticalPeakMass + Properties.peakDifferenceThreshold;
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
				impValue = 1;
				return impValue;
			}
			
			//find out ionMatchTally and intensity distributions
			int totalIonsAbove50 = 0;
			int totalIonsAbove25 = 0;
			int totalIonsAbove12 = 0;
			int totalIonsAbove06 = 0;
			double totalMatchingIntensity = 0.0;
			boolean yIonMatch, bIonMatch, ionMatch = false;
			int acidMatchTally = 0;
			int yGreaterThanBTally = 0;
			for (i = 0; i < peptideLengthMinusOne; i++) {
				yIonMatch = yIonMatchesWithHighestIntensity[i] > 0.0;
				bIonMatch = bIonMatchesWithHighestIntensity[i] > 0.0;
				ionMatch = bIonMatch || yIonMatch;
				if (ionMatch) acidMatchTally++;
				if (yIonMatchesWithHighestIntensity[i] > bIonMatchesWithHighestIntensity[i]) yGreaterThanBTally++;
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
			
			//TODO: optimize this.  it is a great place to improve performance
			if (peakMatchProbability > 0.25) {
				impValue = 1;
				return impValue;
			}
			
			//y greater than b probability
			n = peptideLengthMinusOne;
			k = yGreaterThanBTally;
			double yGreaterThanBProbability = MathFunctions.getCachedBinomialProbability50(n, k);
			
			
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
			
			impValue =  peakMatchProbability  * intensityProbability * yGreaterThanBProbability;
		}

		//this is a normalizing factor as a true match with a long peptide will get a greater
		//score than a true match with a short pepitide, though they are equally true
		impValue *= MathFunctions.cachedLog(peptide.getAcidSequence().length);
		if (impValue > 1) impValue = 1;
		return impValue;
	}
	
	
	public int compareTo(Match match) {
			if (sortParameter == SORT_BY_SCORE) {
				//want to sort from greatest to least
				if (score > match.getScore()) return -1;
				if (score < match.getScore()) return  1;
				return 0;
			} else
			if (sortParameter == SORT_BY_SCORE_RATIO) {
				//want to sort from greatest to least
				if (scoreRatio > match.getScoreRatio()) return -1;
				if (scoreRatio < match.getScoreRatio()) return  1;
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
			if (sortParameter == SORT_BY_E_VALUE) {
				if (eValue < match.getEValue()) return -1;
				if (eValue > match.getEValue()) return  1;
				return 0;
			} else
			if (sortParameter == SORT_BY_IMP_VALUE) {
				if (impValue < match.getImpValue()) return -1;
				if (impValue > match.getImpValue()) return  1;
				return 0;
			} else
			if (sortParameter == SORT_BY_SPECTRUM_ID_THEN_SCORE) {
				if (spectrum.getId() < match.getSpectrum().getId()) return -1;
				if (spectrum.getId() > match.getSpectrum().getId()) return  1;
				//if spectrum is sorted, also sort by tandemFit
				if (score > match.getScore()) return -1;
				if (score < match.getScore()) return  1;
				return 0;
			} else
			if (sortParameter == SORT_BY_SPECTRUM_ID_THEN_PEPTIDE) {
				//first by spectrum ID
				if (spectrum.getId() < match.getSpectrum().getId()) return -1;
				if (spectrum.getId() > match.getSpectrum().getId()) return  1;
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
			if (sortParameter == SORT_BY_RANK_THEN_E_VALUE) {
				if (rank < match.rank) return -1;
				if (rank > match.rank) return  1;
				if (eValue < match.getEValue()) return -1;
				if (eValue > match.getEValue()) return  1;
				return 0;
			} else 
			if (sortParameter == SORT_BY_RANK_THEN_SCORE) {
				if (rank < match.rank) return -1;
				if (rank > match.rank) return  1;
				if (score < match.getScore()) return 1;
				if (score > match.getScore()) return -1;
				return 0;
			} else 
			if (sortParameter == SORT_BY_SPECTRUM_PEPTIDE_MASS_DIFFERENCE) {
				//i'm putting this calculation in here as this is a not-often-used sort
				//so calculating and storing this for every match is unnecessary
				double myDifference = spectrum.getMass() - peptide.getMass();
				double theirDifference = match.getSpectrum().getMass() - match.getPeptide().getMass();
				if (myDifference > theirDifference) return -1;
				if (myDifference < theirDifference) return  1;
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

	public double getImpValue() {
		return impValue;
	}

	public boolean equals(Match match) {
		if (getScore() == match.getScore()) {
			if (peptide.equals(match.getPeptide())) return true;
		}
		return false;
	}
	
	/**
	 * Returns the default score.  This could be TandemFit or HMM.
	 */
	public double getScore() {
		return score;
	}

	
	/**
	 * @return the spectrum
	 */
	public Spectrum getSpectrum() {
		return spectrum;
	}

	public double getScoreRatio() {
		return scoreRatio;
	}

	/**
	 * @return the peptide
	 */
	public Peptide getPeptide() {
		return peptide;
	}

	public double getEValue() {
		return eValue;
	}

	/**
	 * @return the ionMatchTally
	 */
	public int getIonMatchTally() {
		return ionMatchTally;
	}

	public static void setSortParameter(int sortParameter) {
		Match.sortParameter = sortParameter;
	}


	public void setScoreRatio(double scoreRatio) {
		this.scoreRatio = scoreRatio;
	}

	public void setEValue(double eValue) {
		this.eValue = eValue;
	}

	public String toString() {
		return spectrum.getId() + "\t" + peptide.getAcidSequenceString()  + "\t" + getScore() + "\t" + getEValue() + "\t" + getImpValue();
	}
	
	public double calculateEValue() {
		eValue = spectrum.getEValue(getScore());
		return eValue;
	}
	
	public boolean hasModification() {
		return false;
	}
	
	public ArrayList<Modification> getModifications() {
		return new ArrayList<Modification>(peptide.getLength());
	}

}
