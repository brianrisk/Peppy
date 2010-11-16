package Peppy;

import HMMScore.HMMClass;
import Statistics.HasEValue;

/**
 * An object which contains scoring mechanisms to evaluate a spectrum/peptide match.
 * @author Brian Risk
 *
 */
public class Match implements Comparable<Match>, HasEValue{
	
	private double score = 0.0;
	private double scoreTandemFit = 0.0;
	private double tandemFitScoreRatio = 0.0;
	private int rankCount = 0; //longest variable name I've ever made
	private double scoreHMM = 0.0;
	private double eValue;
	public int ionMatchTally = 0;
	private int rank = Integer.MAX_VALUE;
	
	private Spectrum spectrum;
	private Peptide peptide;
	private Sequence sequence;
	
	public final static int SORT_BY_SCORE = 0;
	public final static int SORT_BY_SPECTRUM_ID_THEN_SCORE = 1;
	public final static int SORT_BY_LOCUS = 2;
	public final static int SORT_BY_SCORE_RATIO = 3;
	public final static int SORT_BY_HMM = 4;
	public final static int SORT_BY_TANDEM_FIT = 5;
	public final static int SORT_BY_E_VALUE = 6;
	public final static int SORT_BY_RANK_THEN_E_VALUE = 7;
	public final static int SORT_BY_SPECTRUM_PEPTIDE_MASS_DIFFERENCE = 8;
	public final static int SORT_BY_SPECTRUM_ID_THEN_PEPTIDE = 9;
	public final static int SORT_BY_RANK_THEN_SCORE = 10;
	
	
	
	private static int sortParameter = SORT_BY_SCORE;
	
	public Match(String spectrumString, String peptideString) {
		this.spectrum = new Spectrum(spectrumString);
		this.peptide = new Peptide(peptideString);
		this.sequence = null;
		calculateScore();
	}
	
	public Match(Spectrum spectrum, Peptide peptide, Sequence sequence) {
		this.spectrum = spectrum;
		this.peptide = peptide;
		this.sequence = sequence;
		calculateScore();
	}
	
	public void calculateScore() {
		if (Properties.defaultScore == Properties.DEFAULT_SCORE_TANDEM_FIT) {
			score = calculateTandemFit();
		} else
		if (Properties.defaultScore == Properties.DEFAULT_SCORE_HMM) {
			score = calculateHMM();
		}
	}
	
	public double calculateTandemFit() {
		String peptideString = peptide.getAcidSequence();

		int i;
		boolean atLeastOneMatch = false;
		double theoreticalPeakMass, peakMass;
		int peakIndex, seqIndex;
		
		//we want -1 because most of these spectra will have a match with 
		//the last theoretical peak
		int peptideLengthMinusOne = peptideString.length() - 1;
		
		double [] bIonMatchesWithHighestIntensity = new double[peptideString.length()];
		for (i = 0; i < peptideString.length(); i++) bIonMatchesWithHighestIntensity[i] = 0.0;
		double [] yIonMatchesWithHighestIntensity = new double[peptideString.length()];
		for (i = 0; i < peptideString.length(); i++) yIonMatchesWithHighestIntensity[i] = 0.0;

		
		//find the ranges around our theoretical peptides where we
		//count spectrum peaks
		double [] theoreticalPeaksLeft = new double[peptideLengthMinusOne];
		double [] theoreticalPeaksRight = new double[peptideLengthMinusOne];
		
		
		/* y-ion  */
		//computing the left and right boundaries for the ranges where our peaks should land
		theoreticalPeakMass = peptide.getMass() + Properties.rightIonDifference;
		for (i = 0; i < peptideLengthMinusOne; i++) {
			theoreticalPeakMass -= Definitions.getAminoAcidWeightMono(peptideString.charAt(i));
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
			scoreTandemFit = 0.0;
			return 0.0;
		}
			
		
		/* b-ion  */
		theoreticalPeakMass = Properties.leftIonDifference;
		for (i = 0; i < peptideLengthMinusOne; i++) {
			theoreticalPeakMass += Definitions.getAminoAcidWeightMono(peptideString.charAt(i));
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
			scoreTandemFit = 0.0;
			return 0.0;
		}
		
		//find out final tally
		double score = 0.0;
		boolean yIonTrue, bIonTrue;
		double amountToAdd;
		for (i = 0; i < peptideString.length(); i++) {
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
		scoreTandemFit = score;
		return score;
	}
	
	public double calculateHMM() {
		HMMClass scorer = new HMMClass(peptide.getAcidSequence(), spectrum);
		scoreHMM = scorer.score();
		return scoreHMM;
	}
	
	public int compareTo(Match match) {
			if (sortParameter == SORT_BY_SCORE) {
				//want to sort from greatest to least
				if (score > match.getScore()) return -1;
				if (score < match.getScore()) return  1;
				return 0;
			} else
			if (sortParameter == SORT_BY_HMM) {
				//want to sort from greatest to least
				if (scoreHMM > match.getScoreHMM()) return -1;
				if (scoreHMM < match.getScoreHMM()) return  1;
				return 0;
			} else
			if (sortParameter == SORT_BY_TANDEM_FIT) {
				//want to sort from greatest to least
				if (scoreTandemFit > match.getScoreTandemFit()) return -1;
				if (scoreTandemFit < match.getScoreTandemFit()) return  1;
				return 0;
			} else
			if (sortParameter == SORT_BY_SCORE_RATIO) {
				//want to sort from greatest to least
				if (tandemFitScoreRatio > match.getTandemFitScoreRatio()) return -1;
				if (tandemFitScoreRatio < match.getTandemFitScoreRatio()) return  1;
				return 0;
			} else
			if (sortParameter == SORT_BY_LOCUS) {
				if (sequence.getId() < match.getSequence().getId()) return -1;
				if (sequence.getId() > match.getSequence().getId()) return  1;
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
				return peptide.getAcidSequence().compareTo(match.getPeptide().getAcidSequence());
			} else	
			if (sortParameter == SORT_BY_RANK_THEN_E_VALUE) {
				if (rank < match.getRank()) return -1;
				if (rank > match.getRank()) return  1;
				if (eValue < match.getEValue()) return -1;
				if (eValue > match.getEValue()) return  1;
				return 0;
			} else 
			if (sortParameter == SORT_BY_RANK_THEN_SCORE) {
				if (rank < match.getRank()) return -1;
				if (rank > match.getRank()) return  1;
				if (score < match.getScore()) return 1;
				if (score > match.getScore()) return -1;
				return 0;
			} else 
			if (sortParameter == SORT_BY_SPECTRUM_PEPTIDE_MASS_DIFFERENCE) {
				//i'm putting this calculation in here as this is a not-often-used sort
				//so calculating and storing this for every match is unnecessary
				double myDifference = spectrum.getPrecursorMass() - peptide.getMass();
				double theirDifference = match.getSpectrum().getPrecursorMass() - match.getPeptide().getMass();
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

	public boolean equals(Match match) {
		if (getScore() == match.getScore()) {
			if (peptide.equals(match.getPeptide())) return true;
		}
		return false;
	}

	public int getRank() {
		return rank;
	}

	public void setRank(int rank) {
		this.rank = rank;
	}

	public int getRankCount() {
		return rankCount;
	}

	public void setRankCount(
			int rankCount) {
		this.rankCount = rankCount;
	}

	public double getScoreTandemFit() {
		return scoreTandemFit;
	}
	
	
	/**
	 * Returns the default score.  This could be TandemFit or HMM.
	 */
	public double getScore() {
		if (Properties.defaultScore == Properties.DEFAULT_SCORE_TANDEM_FIT) {
			return scoreTandemFit;
		}
		return scoreHMM;
	}

	public double getScoreHMM() {
		return scoreHMM;
	}
	
	/**
	 * @return the spectrum
	 */
	public Spectrum getSpectrum() {
		return spectrum;
	}

	public double getTandemFitScoreRatio() {
		return tandemFitScoreRatio;
	}

	/**
	 * @return the peptide
	 */
	public Peptide getPeptide() {
		return peptide;
	}
	
	public Sequence getSequence() {
		return sequence;
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

	/**
	 * @param spectrum the spectrum to set
	 */
	public void setSpectrum(Spectrum spectrum) {
		this.spectrum = spectrum;
	}

	public void setTandemFitScoreRatio(double tandemFitScoreRatio) {
		this.tandemFitScoreRatio = tandemFitScoreRatio;
	}

	/**
	 * @param peptide the peptide to set
	 */
	public void setPeptide(Peptide peptide) {
		this.peptide = peptide;
	}

	public void setEValue(double eValue) {
		this.eValue = eValue;
	}

	public String toString() {
		return spectrum.getId() + "\t" + peptide.getAcidSequence()  + "\t" + getScore();
	}
	
	public double calculateEValue() {
		eValue = spectrum.getEValue(getScore());
		return eValue;
	}

}
