package Statistics;

import java.util.ArrayList;

import Utilities.U;

/**
 * For the caluclation and storing of E Values.
 * 
 * This will do all of the housekeeping things like
 * keep track of histograms, find least squares fit, etc.
 * 
 * @author Brian Risk
 *
 */

public class EValueCalculator {
	
	//E-Values:  allocating histogram variables
	private final int numberOfHistogramBars = 100;
	private int [] histogram = new int[numberOfHistogramBars];
	private double [] scoreProbabilities = new double[numberOfHistogramBars];
	private double [] survivability = new double[numberOfHistogramBars];
	private double [] xValues = new double[numberOfHistogramBars];
	private double highScore = -1.0;
	private double lowScore = 0.0;
	private double barWidth;
	private double m;
	private double b;
	private int numberOfMatches = 0;
	
	public EValueCalculator(ArrayList<? extends HasEValue> values, ArrayList<? extends HasEValue> topValues) {
		addScores(values, topValues);
	}
	
	public EValueCalculator(int [] histogram) {
		this.histogram = histogram;
	}
	
	/**
	 * find expected value (a.k.a. "e value") for top matches
	 * 
	 * This method assumes that matchesForOneSpectrum is already sorted from highest score to lowest.
	 * Calculates e values for each of the top matches
	 * @param matchesForOneSpectrum
	 * @param topMatches
	 */
	public void addScores(ArrayList<? extends HasEValue> values, ArrayList<? extends HasEValue> topValues) {
		//Set these values if this is our first time calculating e value
		if (highScore < 0) {
			//setting up the histogram parameters
			highScore = topValues.get(0).getScore();
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
		numberOfMatches += values.size();

		//populate the histogram
		int bin;
		for (HasEValue value: values) {
			bin = (int) Math.floor((value.getScore() - lowScore) / barWidth);
			if (bin < numberOfHistogramBars) {
				histogram[bin]++;
			} else {
				histogram[numberOfHistogramBars - 1]++;
			}
		}
		
		//find score probabilities
		for (int i = 0; i < numberOfHistogramBars; i++) {
			scoreProbabilities[i] = (double) histogram[i] / numberOfMatches;
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
		chopIndex = 0;
		topIndex = numberOfHistogramBars;
		
		//taking the log of each of the survivability.  Only concerned
		//about values at and above chopIndex
		for (int i = chopIndex; i < topIndex; i++) {
			survivability[i] = Math.log(survivability[i]);
		}
		
		//finding the least squares fit for that region
		// y = m * x + b
		m = U.calculateM(xValues, survivability, histogram, chopIndex, topIndex);
		b = U.calculateB(xValues, survivability, histogram, chopIndex, topIndex, m);
		
		//using our m and b to derive e values for all top matches
		double eValue;
		for (HasEValue value: topValues) {
			eValue = getEValue(value.getScore());
			value.setEValue(eValue);
		}
	}

	/**
	 * find the e value for just one given score.
	 * Returns largest double value if result is NaN
	 * @param score
	 * @return
	 */
	public double getEValue(double score) {
		double eValue = m * score + b;
		eValue = Math.exp(eValue);
		eValue *= numberOfMatches;
		//setting to Double.MAX_VALUE if eValue is Nan
		if (eValue <= 1 == eValue >= 1) eValue = Double.MAX_VALUE;
		return eValue;
//		return getNormalDistributionEValue(score);
	}


	public int [] getHistogram() {return histogram;}
	
	public void setHistogram(int[] histogram) {
		this.histogram = histogram;
	}
	
	/* 
	 * E value experiments
	 */
	double mean = -1;
	double variance = -1;
	int numberOfHits = -1;
	double scalar = -1;
	private void calculateDistribution() {
		//calculate mean
		double total = 0;
		numberOfHits = 0;
		for (int i = 0; i < numberOfHistogramBars; i++) {
			numberOfHits += histogram[i];
			total += xValues[i] * histogram[i];
		}
		mean = total / numberOfHits;
		
		//calculate variance
		double difference;
		variance = 0;
		for (int i = 0; i < numberOfHistogramBars; i++) {
			difference = xValues[i] - mean;
			variance += histogram[i] * difference * difference;
		}
		variance /= numberOfHits;
		
		scalar = 1.0 / Math.sqrt(2.0 * Math.PI * variance);
		
	}
	
	public double getNormalDistributionEValue(double score) {
		calculateDistribution();
		
		//exponent
		double numerator = (score - mean);
		numerator *= numerator;
		double exponent =  numerator / (-2.0 * variance);
		
		return (double)  numberOfHits * scalar * Math.exp(exponent);
	}
	

}
