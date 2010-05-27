package Peppy;

import HMMScore.HMMClass;
	import Utilities.U;

public class SpectrumPeptideMatch implements Comparable<SpectrumPeptideMatch>{
	
	private double score = 0.0;
	private double scoreMSMSFit = 0.0;
	private double MSMSFitScoreRatio = 0.0;
	private int MSMSFitRank = -1;
	private double scoreHMM = 0.0;
	private double eValue;
	
	private Spectrum spectrum;
	private Peptide peptide;
	private Sequence sequence;
	
	final static double useAcidThreshold = 100.0;
	
	public final static int DEFAULT_SCORE_MSMS_FIT = 0;
	public final static int DEFAULT_SCORE_HMM = 1;
	private static int defaultScore = DEFAULT_SCORE_MSMS_FIT;
//	private static int defaultScore = DEFAULT_SCORE_HMM;
	
	public final static int SORT_BY_DEFAULT = 0;
	public final static int SORT_BY_SPECTRUM_ID = 1;
	public final static int SORT_BY_LOCUS = 2;
	public final static int SORT_BY_SCORE_RATIO = 3;
	public final static int SORT_BY_HMM = 4;
	public final static int SORT_BY_MSMSFIT = 5;
	public final static int SORT_BY_E_VALUE = 6;
	private static int sortParameter = SORT_BY_DEFAULT;
	
	public SpectrumPeptideMatch(Spectrum spectrum, Peptide peptide, Sequence sequence) {
		this.spectrum = spectrum;
		this.peptide = peptide;
		this.sequence = sequence;
		calculateScore();
	}
	
	public void calculateScore() {
//		calculateMSMSFit();
		if (defaultScore == DEFAULT_SCORE_MSMS_FIT) {
			score = calculateMSMSFit();
		} else
		if (defaultScore == DEFAULT_SCORE_HMM) {
			score = calculateHMM();
		}
	}
	
	public double calculateMSMSFit() {
		String peptideString = peptide.getSequence();

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
			scoreMSMSFit = 0.0;
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
			scoreMSMSFit = 0.0;
			return 0.0;
		}
		
		//find out final tally
		double score = 0.0;
		for (i = 0; i < peptideString.length(); i++) {
			if (yIonMatchesWithHighestIntensity[i] > 0.0) score += Math.pow(yIonMatchesWithHighestIntensity[i], Properties.peakIntensityExponent);
			if (bIonMatchesWithHighestIntensity[i] > 0.0) score += Math.pow(bIonMatchesWithHighestIntensity[i], Properties.peakIntensityExponent);
		}
		scoreMSMSFit = score;
		return score;
	}
	
	public double calculateHMM() {
		HMMClass scorer = new HMMClass(peptide.getSequence(), spectrum);
		scoreHMM = scorer.score();
		return scoreHMM;
	}

	public int compareTo(SpectrumPeptideMatch match) {
		if (sortParameter == SORT_BY_HMM) {
			//want to sort from greatest to least
			if (scoreHMM > match.getScoreHMM()) return -1;
			if (scoreHMM < match.getScoreHMM()) return  1;
			return 0;
		} else
		if (sortParameter == SORT_BY_MSMSFIT) {
			//want to sort from greatest to least
			if (scoreMSMSFit > match.getScoreMSMSFit()) return -1;
			if (scoreMSMSFit < match.getScoreMSMSFit()) return  1;
			return 0;
		} else
		if (sortParameter == SORT_BY_SCORE_RATIO) {
			//want to sort from greatest to least
			if (MSMSFitScoreRatio > match.getMSMSFitScoreRatio()) return -1;
			if (MSMSFitScoreRatio < match.getMSMSFitScoreRatio()) return  1;
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
			//if spectrum is sorted, also sort by msmsfit
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

	public double getScoreMSMSFit() {
		return scoreMSMSFit;
	}
	
	
	/**
	 * Returns the default score.  This could be MSMSFit or HMM.
	 */
	public double getScore() {
		if (defaultScore == DEFAULT_SCORE_MSMS_FIT) {
			return scoreMSMSFit;
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

	public double getMSMSFitScoreRatio() {
		return MSMSFitScoreRatio;
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

	public int getMSMSFitRank() {
		return MSMSFitRank;
	}

	public double getEValue() {
		return eValue;
	}

	public static void setSortParameter(int sortParameter) {
		SpectrumPeptideMatch.sortParameter = sortParameter;
	}

	/**
	 * @param spectrum the spectrum to set
	 */
	public void setSpectrum(Spectrum spectrum) {
		this.spectrum = spectrum;
	}

	public void setMSMSFitScoreRatio(double mSMSFitScoreRatio) {
		MSMSFitScoreRatio = mSMSFitScoreRatio;
	}

	/**
	 * @param peptide the peptide to set
	 */
	public void setPeptide(Peptide peptide) {
		this.peptide = peptide;
	}

	
	public void setMSMSFitRank(int mSMSFitRank) {
		MSMSFitRank = mSMSFitRank;
	}

	public void setEValue(double eValue) {
		this.eValue = eValue;
	}

	public String toString() {
		return peptide.getSequence() + " " + scoreMSMSFit + " " + peptide.getMass();
	}

}
