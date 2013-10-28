package Peppy;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

import Graphs.HistogramVisualizer;
import Math.LeastSquares;

/**
 * We want to hold collections of matches for various things.  All matches for a given peptide,
 * or spectrum or region or protein... This abstract class defines common functionality for all
 * of these.
 * 
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public abstract class Matches implements Comparable<Matches> {
	
	ArrayList<Match> matches = new ArrayList<Match>();
	
	/* the top score of all of our matches */
	double score = Properties.minimumScore;
		
	
	private static int labelTracker = 0;
	public static final int KEEP_ONLY_BEST_MATCHES = labelTracker++;
	public static final int KEEP_MATCHES_AT_MINIMUM_SCORE = labelTracker++;
	
	private int whatToKeep = KEEP_ONLY_BEST_MATCHES;
	
	/* only keep decoy hits if they trump our existing score 
	 * For FDR we only need to know when the decoy trumps the target hits */
	private boolean ignoreLesserDecoys = true;
	
	/* properties for calculating e values */
	private final int numberOfHistogramBars = 1;
	private int [] histogram = new int[numberOfHistogramBars];
	private int [] smoothedHistogram  = new int[numberOfHistogramBars];
	private double [] scoreProbabilities = new double[numberOfHistogramBars];
	private double [] survivability = new double[numberOfHistogramBars];
	private double m;
	private double b;
	private int numberOfMatches = 0;
	
	/**
	 * Add a match only if it's score is greater than or equal to the reigning score.
	 * If it is greater than, then the existing matches are cleared out.
	 * 
	 * NOTE: don't add if it is a duplicate match??
	 * 
	 * @param match
	 */
	public void addMatch(Match match) {
		/* tracking for eValues */
//		numberOfMatches++;
//		int histogramIndex = (int) Math.round(match.getScore() * 2);
//		if (histogramIndex >= numberOfHistogramBars) histogramIndex = numberOfHistogramBars -1;
//		U.p(histogramIndex);
//		histogram[histogramIndex]++;
		
		
		/* selecting which matches should be added */
		if (ignoreLesserDecoys) {
			if (match.getPeptide().isDecoy()) {
				if (match.getScore() <= score) {
					return;
				} 
			}
		}
		
		if (whatToKeep == KEEP_ONLY_BEST_MATCHES) {
			if (match.getScore() > score) {
				matches.clear();
				matches.add(match);
				score = match.getScore();
			} else {
				if (match.getScore() == score) {
					
					/* For DNA there may be many, many matches.  This is saying, that, if there
					 * are more than four matches, then the efficacy of this result is low.
					 * Therefore, the memory and time to track each match is not necessary.
					 * 
					 * This one line of code can save huge amounts of memory.
					 */
					if (Properties.isSequenceFileDNA && matches.size() >= 4) return;
					
					matches.add(match);
				}
			}
			return;
		}
		

		if (whatToKeep == KEEP_MATCHES_AT_MINIMUM_SCORE) {
			if (match.getScore() >= Properties.minimumScore) {
				matches.add(match);
				if (match.getScore() > score) {
					score = match.getScore();
				}
			} 
		}
	}
	
	public void clearMatches() {
		matches.clear();
		score = Properties.minimumScore;
	}
	
	public ArrayList<Match> getMatches() {
		return matches;
	}

	public double getScore() {
		if (matches.size() > 0) {
			return score;
		} else {
			return 0;
		}
	}
	
	public static ArrayList<Match> getBestMatches(ArrayList<Match> matches) {
		ArrayList<Match> bestMatches = new ArrayList<Match>();
		Match.setSortParameter(Match.SORT_BY_PEPTIDE_THEN_SCORE);
		Collections.sort(matches);
		Peptide peptide = new Peptide("k");
		for (Match match: matches) {
			if (!peptide.equals(match.getPeptide())) {
				peptide = match.getPeptide();
				bestMatches.add(match);
			}
		}
		Match.setSortParameter(Match.SORT_BY_SCORE);
		Collections.sort(matches);
		Collections.sort(bestMatches);
		return bestMatches;
	}
	
	public static ArrayList<Match> getMatchesWithSpectrum(Spectrum spectrum, ArrayList<Match> theseMatches) {
		ArrayList<Match> out = new ArrayList<Match>();
		for (int i = 0; i < theseMatches.size(); i++) {
			Match match = theseMatches.get(i);
			if (match.getSpectrum().getId() == spectrum.getId()) {
				out.add(match);
			}
		}
		return out;
	}
	
	public static ArrayList<Match> getMatchesWithSequence(Sequence_DNA sequence_DNA, ArrayList<Match> theseMatches) {
		ArrayList<Match> out = new ArrayList<Match>();
		for (int i = 0; i < theseMatches.size(); i++) {
			Match match = theseMatches.get(i);
			if (match.getPeptide().getParentSequence().getId() == sequence_DNA.getId()) {
				out.add(match);
			}
		}
		return out;
	}

	public int getWhatToKeep() {
		return whatToKeep;
	}

	public void setWhatToKeep(int whatToKeep) {
		this.whatToKeep = whatToKeep;
	}
	
	
	public int compareTo(Matches o) {
		if (getScore() < o.getScore()) return 1;
		if (getScore() > o.getScore()) return -1;
		return 0;
	}

	public boolean ignoreLesserDecoys() {
		return ignoreLesserDecoys;
	}

	public void setIgnoreLesserDecoys(boolean ignoreLesserDecoys) {
		this.ignoreLesserDecoys = ignoreLesserDecoys;
	}
	
	
	private void calculateHistogramProperties() {
		smoothedHistogram = smoothHistogram(histogram, 1);
//		smoothedHistogram = histogram;
		
		//find score probabilities
		for (int i = 0; i < numberOfHistogramBars; i++) {
			scoreProbabilities[i] = (double) smoothedHistogram[i] / numberOfMatches;
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
		for (topIndex = chopIndex + 1; topIndex < numberOfHistogramBars; topIndex++) {
			if (smoothedHistogram[topIndex] == 0) break;
		}
		if (topIndex >= numberOfHistogramBars) topIndex = numberOfHistogramBars;
		
		//if width is only one bar wide, make it 2
		if (topIndex - chopIndex == 1) chopIndex--;
		if (chopIndex < 0) chopIndex = 0;
		
		//taking the log of each of the survivability.  Only concerned
		//about values at and above chopIndex
		for (int i = chopIndex; i < topIndex; i++) {
			survivability[i] = Math.log(survivability[i]);
		}
		
		//finding the least squares fit for that region
		// y = m * x + b
		m = LeastSquares.calculateM(survivability, smoothedHistogram, chopIndex, topIndex);
		b = LeastSquares.calculateB( survivability, smoothedHistogram, chopIndex, topIndex, m);
		
		
//		try {
//			U.p(m + ", " + b);
//			U.p("drawing...");
//			HistogramVisualizer.drawHistogram(histogram, 500, 1500, new File("histogram.jpg"));
//		} catch (IOException e) {
//			U.p("problem happened");
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
	}
	
	
	public int[] smoothHistogram(int [] histogram, int radius) {
		int [] out = new int[histogram.length];
		double sigma = radius / 3.0;
		double coefficient = 1.0 / Math.sqrt(2 * Math.PI * sigma * sigma);
		double denominator = -2 * sigma * sigma;
		int gaussianSize = radius * 2 + 1;
		double [] gaussianY = new double[gaussianSize];
		int [] gaussianX = new int[gaussianSize];
		for (int i = 0; i < gaussianSize; i++) {
			gaussianX[i] = i - radius;
			gaussianY[i] = coefficient * Math.exp(gaussianX[i] * gaussianX[i] / denominator);
		}
		double smoothedValue;
		int index;
		for (int i = 0; i < histogram.length; i++) {
			smoothedValue = 0;
			for (int j = 0; j < gaussianSize; j++) {
				index = gaussianX[j] + i;
				if (index < 0) continue;
				if (index >= histogram.length) break;
				smoothedValue += gaussianY[j] * histogram[index];
			}
			out[i] = (int) Math.round(smoothedValue);
		}
		return out;
	}

	/**
	 * find the e value for just one given score.
	 * Returns largest double value if result is NaN
	 * @param score
	 * @return
	 */
	private double calculateEValueOfScore(double score) {
		double eValue = m * score + b;
		eValue = Math.exp(eValue);
		eValue *= numberOfMatches;
		//setting to Double.MAX_VALUE if eValue is Nan
		if (eValue <= 1 == eValue >= 1) eValue = Double.MAX_VALUE;
		return eValue;
	}
	
	public void calculateEValues() {
		calculateHistogramProperties();
		score = Double.MIN_VALUE;
		double eValue;
		for (Match match: matches) {
			eValue = calculateEValueOfScore(match.getScore());
			eValue = -1.0 * Math.log10(eValue);
//			if(!match.getPeptide().isDecoy()) U.p(eValue);
			match.setScore(eValue);
			if (score < eValue) score = eValue;
		}
	}

}
