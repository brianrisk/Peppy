package Peppy;
import java.util.Collections;
import java.util.ArrayList;
import java.util.Vector;

import Utilities.U;


public class ScoringThread implements Runnable {
	
	ArrayList<Peptide> peptides;
	Spectrum spectrum;
	ScoringEngine scoringEngine;
	Sequence sequence;
	
	//E-Values:  allocating histogram variables
	int numberOfHistogramBars = 50;
	int [] histogram = new int[numberOfHistogramBars];
	double [] scoreProbabilities = new double[numberOfHistogramBars];
	double [] survivability = new double[numberOfHistogramBars];
	double [] xValues = new double[numberOfHistogramBars];
	
	/**
	 * @param peptides
	 * @param spectrum
	 */
	public ScoringThread(Spectrum spectrum, ArrayList<Peptide> peptides, ScoringEngine scoringEngine, Sequence sequence) {
		this.spectrum = spectrum;
		this.peptides = peptides;
		this.scoringEngine = scoringEngine;
		this.sequence = sequence;
	}

	public void run() {
		
		while (spectrum != null) {
			
			ArrayList<SpectrumPeptideMatch> matchesForOneSpectrum = new ArrayList<SpectrumPeptideMatch>();
	
			//find the first index of the peptide with mass greater than lowestPeptideMassToConsider
			double lowestPeptideMassToConsider = spectrum.getPrecursorMass() - Properties.spectrumToPeptideMassError;
			int firstPeptideIndex = findFirstIndexWithGreaterMass(peptides, lowestPeptideMassToConsider);
			
			//find the first index of the peptide with mass greater than highestPeptideMassToConsider
			double highestPeptideMassToConsider = spectrum.getPrecursorMass() + Properties.spectrumToPeptideMassError;
			int lastPeptideIndex = findFirstIndexWithGreaterMass(peptides, highestPeptideMassToConsider);
			
			//examine only peptides in our designated mass range
			for (int peptideIndex = firstPeptideIndex; peptideIndex < lastPeptideIndex; peptideIndex++) {
				Peptide peptide = peptides.get(peptideIndex);
				SpectrumPeptideMatch match = new SpectrumPeptideMatch(spectrum, peptide, sequence);
				if (match.getScore() == 0.0) {
					continue;
				}
				matchesForOneSpectrum.add(match);
			}
			
			//collect the top maximumNumberOfMatchesForASpectrum
			Collections.sort(matchesForOneSpectrum);
			ArrayList<SpectrumPeptideMatch> topMatches = new ArrayList<SpectrumPeptideMatch>();
			int max = Properties.maximumNumberOfMatchesForASpectrum;
			if (matchesForOneSpectrum.size() < max) max = matchesForOneSpectrum.size();
			for (int i = 0; i < max; i++) {
				topMatches.add(matchesForOneSpectrum.get(i));
			}
			
			//assign E values to top Matches:
			if (matchesForOneSpectrum.size() == 0) {
				//U.p("There were zero matches for this spectrum file: " + spectrum.getFile().getName());
			} else {
				calculateEValues(matchesForOneSpectrum, topMatches);
			}
			
			//return results, get new task
			spectrum = scoringEngine.getNextSpectrum(topMatches);
		}
	}
	
	/**
	 * This method assumes that matchesForOneSpectrum is already sorted from highest score to lowest.
	 * @param matchesForOneSpectrum
	 * @param topMatches
	 */
	private void calculateEValues(ArrayList<SpectrumPeptideMatch> matchesForOneSpectrum, ArrayList<SpectrumPeptideMatch> topMatches) {
		/*
		 * find expected value (aka "e value") for top matches
		 */
		//setting up the histogram paramaters
		double highScore = matchesForOneSpectrum.get(0).getScore();
//		double lowScore = matchesForOneSpectrum.get(matchesForOneSpectrum.size() - 1).getScoreMSMSFit();
		double lowScore = 0;
		double barWidth = (highScore - lowScore) / numberOfHistogramBars;
		
		//initializing histograms and xValues
		for (int i = 0; i < numberOfHistogramBars; i++) {
			histogram[i] = 0;
			xValues[i] = lowScore + (i * barWidth);
		}

		//populate the histogram
		int bin;
		for (SpectrumPeptideMatch match: matchesForOneSpectrum) {
			bin = (int) Math.floor((match.getScore() - lowScore) / barWidth);
			if (bin == numberOfHistogramBars) bin = numberOfHistogramBars - 1;
			histogram[bin]++;
		}
		
		//find score probabilities
		for (int i = 0; i < numberOfHistogramBars; i++) {
			scoreProbabilities[i] = (double) histogram[i] / matchesForOneSpectrum.size();
		}
		
		//find survivability values
		survivability[numberOfHistogramBars - 1] = scoreProbabilities[numberOfHistogramBars - 1];
		for (int i = numberOfHistogramBars - 2; i >= 0; i--) {
			survivability[i] = survivability[i + 1] + scoreProbabilities[i];
		}
		
		//find index survivability values at 0.1 or less
		int chopIndex;
		for (chopIndex = 0; chopIndex < numberOfHistogramBars; chopIndex++) {
			if (survivability[chopIndex] <= 0.1) break;
		}
		
		//find first 0 above chopIndex
		int topIndex;
		for (topIndex = chopIndex; topIndex < numberOfHistogramBars; topIndex++) {
			if (histogram[topIndex] == 0) break;
		}
		//if we don't want to use topIndex....
//		topIndex = numberOfHistogramBars;
		
		//taking the log of each of the survivability.  Only concerned
		//about values at and above chopIndex
		for (int i = chopIndex; i < topIndex; i++) {
			survivability[i] = Math.log(survivability[i]);
		}
		
		//finding the least squares fit for that region
		// y = m * x + b
		double m = calculateM(xValues, survivability, chopIndex, topIndex);
		double b = calculateB(xValues, survivability, chopIndex, topIndex, m);
		
		//using our m and be to derive e values for all top matches
		double eValue;
		int peptideCount = matchesForOneSpectrum.size();
		for (SpectrumPeptideMatch match: topMatches) {
			eValue = m * match.getScore() + b;
			eValue = Math.exp(eValue);
			eValue *= peptideCount;
			match.setEValue(eValue);
		}
	}
	
	private double calculateM(double [] xValues, double [] yValues, int start, int stop) {
		double numerator1, numerator2, denomenator1, denomenator2;
		double numerator = 0.0, denomenator = 0.0;
		double temp1 = 0.0, temp2 = 0.0, temp = 0.0;
		double parameterM;
		int i;
		for (i = start; i < stop; i++) {
			temp += (xValues[i] * yValues[i]);
		}
		numerator1 = (stop - start) * (temp);
		for (i = start; i < stop; i++) {
			temp1 += xValues[i];
		}
		for (i = start; i < stop; i++) {
			temp2 += yValues[i];
		}
		numerator2 = temp1 * temp2;
		numerator = numerator1 - numerator2;
		temp1 = temp2 = 0.0;
		for (i = start; i < stop; i++) {
			temp1 = xValues[i];
			temp2 += (temp1 * temp1);
		}
		denomenator1 = (stop - start) * temp2;

		temp1 = 0.0; 
		for (i = start; i < stop; i++) {
			temp1 += xValues[i];
		}
		denomenator2 = (temp1 * temp1);
		denomenator = denomenator1 - denomenator2;
		parameterM = numerator / denomenator;
		return parameterM;
	}
	
	private double calculateB(double [] xValues, double [] yValues, int start, int stop, double m) {
		double parameterB;
		double temp1, temp2;
		int i;

		temp1 = temp2 = 0.0;
		for (i = start; i < stop; i++) {
			temp1 += xValues[i];
		}

		for (i = start; i < stop; i++) {
			temp2 += yValues[i];
		}
		parameterB = (1.0 / (stop - start)) * (temp2 - m * temp1);
		return parameterB;
	}
	
	/**
	 * Boolean search to locate the first peptide in the SORTED list of peptides that has
	 * a mass greater than the "mass" parameter.
	 * @param peptides
	 * @param mass
	 * @return
	 */
	public static int findFirstIndexWithGreaterMass(ArrayList<Peptide> peptides, double mass) {
		Peptide peptide;
		int index = peptides.size() / 2;
		int increment = index / 2;
		while (increment > 0) {
			peptide = peptides.get(index);
			if (peptide.getMass() > mass) {index -= increment;}
			else {index += increment;}
			increment /= 2;
		}
		return index;
	}
	
	public static int findFirstIndexWithGreaterMass(Peptide [] peptides, double mass) {
		Peptide peptide;
		int index = peptides.length / 2;
		int increment = index / 2;
		while (increment > 0) {
			peptide = peptides[index];
			if (peptide.getMass() > mass) {index -= increment;}
			else {index += increment;}
			increment /= 2;
		}
		return index;
	}

}
