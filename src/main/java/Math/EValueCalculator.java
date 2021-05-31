package Math;

import Peppy.Match;

import java.util.ArrayList;

/**
 * For the calculation and storing of E Values.
 * 
 * This will do all of the housekeeping things like
 * keep track of histograms, find least squares fit, etc.
 * 
 * Copyright 2013, Brian Risk
 * 
 * @author Brian Risk
 *
 */

public class EValueCalculator {
	
	//E-Values:  allocating histogram variables
	private final int numberOfHistogramBars = 100;
	private int [] histogram = new int[numberOfHistogramBars];
	private int [] smoothedHistogram  = new int[numberOfHistogramBars];
	private double [] scoreProbabilities = new double[numberOfHistogramBars];
	private double [] survivability = new double[numberOfHistogramBars];
	private double [] xValues = new double[numberOfHistogramBars];
	private double highScore = -1.0;
	private double lowScore = 0.0;
	private double barWidth;
	private double m;
	private double b;
	private int numberOfMatches = 0;
	
	public int getNumberOfMatches() {
		return numberOfMatches;
	}
	
	public EValueCalculator() {
		//nothing
	}
	
	
	/**
	 * find expected value (a.k.a. "e value") for top matches
	 * 
	 * This method assumes that matchesForOneSpectrum is already sorted from highest score to lowest.
	 * Calculates e values for each of the top matches
	 * @param matchesForOneSpectrum
	 * @param topMatches
	 */
	public void addScores(ArrayList<Match> matches) {
		//Set these values if this is our first time calculating e value
		if (highScore < 0) {
			
			/* Find best score */
			highScore = 0;
			for (Match match: matches) {
				if (match.getScore() > highScore) {
					highScore = match.getScore();
				}
			}
			
			//multiplying high score by 2 as there may be higher scores in other chromosomes
			highScore *= 2;
			barWidth = (highScore - lowScore) / numberOfHistogramBars;
			
			//initializing histograms and xValues
			for (int i = 0; i < numberOfHistogramBars; i++) {
				histogram[i] = 0;
				xValues[i] = lowScore + (i * barWidth);
			}
		}
		
		//add to our tally of matches we've observed
		numberOfMatches += matches.size();

		//populate the histogram
		int bin;
		for (Match match: matches) {
			bin = (int) Math.floor((match.getScore() - lowScore) / barWidth);
			if (bin < numberOfHistogramBars) {
				histogram[bin]++;
			} else {
				histogram[numberOfHistogramBars - 1]++;
			}
		}
		
	}
	
	public void calculateHistogramProperties() {
//		smoothedHistogram = smoothHistogram(histogram, 2);
		smoothedHistogram = histogram;
		
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
		m = LeastSquares.calculateM( survivability, smoothedHistogram, chopIndex, topIndex);
		b = LeastSquares.calculateB( survivability, smoothedHistogram, chopIndex, topIndex, m);
		
	}

	/**
	 * find the e value for just one given score.
	 * Returns largest double value if result is NaN
	 * @param score
	 * @return
	 */
	public double calculateEValueOfScore(double score) {
		double eValue = m * score + b;
		eValue = Math.exp(eValue);
		eValue *= numberOfMatches;
		//setting to Double.MAX_VALUE if eValue is Nan
		if (eValue <= 1 == eValue >= 1) eValue = Double.MAX_VALUE;
		return eValue;
	}

	public int [] getSmoothedHistogram() {return smoothedHistogram;}
	
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
	

}
