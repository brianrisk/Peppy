package Peppy;

import HMMScore.HMMClass;

/**
 * An object which contains scoring mechanisms to evaluate a spectrum/peptide match.
 * @author Brian Risk
 *
 */
public class Match implements Comparable<Match>{
	
	private double score = 0.0;
	private double scoreTandemFit = 0.0;
	private double tandemFitScoreRatio = 0.0;
	private int tandemFitRank = -1;
	private double scoreHMM = 0.0;
	private double eValue;
	public int ionMatchTally = 0;
	
	private Spectrum spectrum;
	private Peptide peptide;
	private Sequence sequence;
	
	final static double useAcidThreshold = 100.0;
	
	public final static int SORT_BY_DEFAULT = 0;
	public final static int SORT_BY_SPECTRUM_ID = 1;
	public final static int SORT_BY_LOCUS = 2;
	public final static int SORT_BY_SCORE_RATIO = 3;
	public final static int SORT_BY_HMM = 4;
	public final static int SORT_BY_TANDEM_FIT = 5;
	public final static int SORT_BY_E_VALUE = 6;
	private static int sortParameter = SORT_BY_DEFAULT;
	
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
		theoreticalPeakMass = peptide.getMass() + Properties.yIonDifference;
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
		theoreticalPeakMass = Properties.bIonDifference;
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
		boolean yIonTrue;
		double amountToAdd;
		for (i = 0; i < peptideString.length(); i++) {
			yIonTrue = yIonMatchesWithHighestIntensity[i] > 0.0;
			if (yIonTrue) {
				score += Math.pow(yIonMatchesWithHighestIntensity[i], Properties.peakIntensityExponent);
				ionMatchTally++;
			}
			if (bIonMatchesWithHighestIntensity[i] > 0.0) {
				//Let's try scaling the b ion weight by a bit
				amountToAdd =  0.9 * Math.pow(bIonMatchesWithHighestIntensity[i], Properties.peakIntensityExponent);
				//want to diminish this if this peak is already a yIon
				if (yIonTrue) amountToAdd *= 1.2;
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
			if(peptide.getIndex() < match.getPeptide().getIndex()) return -1;
			if(peptide.getIndex() > match.getPeptide().getIndex()) return  	1;
			return 0;
		} else
		if (sortParameter == SORT_BY_E_VALUE) {
			if (eValue < match.getEValue()) return -1;
			if (eValue > match.getEValue()) return  1;
			return 0;
		} else
		if (sortParameter == SORT_BY_SPECTRUM_ID) {
			if (spectrum.getId() < match.getSpectrum().getId()) return -1;
			if (spectrum.getId() > match.getSpectrum().getId()) return  1;
			//if spectrum is sorted, also sort by tandemFit
			if (score > match.getScore()) return -1;
			if (score < match.getScore()) return  1;
			return 0;
		} else {
			//we want to sort from greatest to least great
			//so -1 is returned where 1 usually is
			if (score > match.getScore()) return -1;
			if (score < match.getScore()) return  1;
			return 0;
		}
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

	public int getTandemFitRank() {
		return tandemFitRank;
	}

	public double getEValue() {
		return eValue;
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

	
	public void setTandemFitRank(int tandemFitRank) {
		this.tandemFitRank = tandemFitRank;
	}

	public void setEValue(double eValue) {
		this.eValue = eValue;
	}

	public String toString() {
		return peptide.getAcidSequence() + " " + scoreTandemFit + " " + peptide.getMass();
	}

}
