//package Statistics;
//
//import java.util.ArrayList;
//
//import Peppy.Match;
//import Utilities.U;
//
//public class StandardDeviation {
//	
//	private int [] histogram = new int[numberOfHistogramBars];
//	
//	//standard deviation
//	private double score_total = -1.0;
//	private double score_mean = -1.0;
//	private double score_variance = -1.0;
//	private double score_standard_deviation = -1.0;
//	
//	public void calculateEValues(ArrayList<Match> matchesForOneSpectrum, ArrayList<Match> topMatches) {
//		//Set these values if this is our first time calculating e value
//		if (highScore < 0) {
//			//setting up the histogram parameters
//			highScore = matchesForOneSpectrum.get(0).getScore();
//			//multiplying high score by 2 as there may be higher scores in other chromsomes
//			highScore *= 2;
//			barWidth = (highScore - lowScore) / numberOfHistogramBars;
//			
//			//initializing histograms and xValues
//			for (int i = 0; i < numberOfHistogramBars; i++) {
//				histogram[i] = 0;
//				xValues[i] = lowScore + (i * barWidth);
//			}
//		}
//		
//		//add to our tally of matches we've observed
//		numberOfMatches += matchesForOneSpectrum.size();
//
//		//populate the histogram
//		int bin;
//		for (Match match: matchesForOneSpectrum) {
//			bin = (int) Math.floor((match.getScore() - lowScore) / barWidth);
//			if (bin < numberOfHistogramBars) {
//				histogram[bin]++;
//			} else {
//				histogram[numberOfHistogramBars - 1]++;
//			}
//		}
//
//		//ha ha, jump in here and find standard deviation distance
//		peptideCount += matchesForOneSpectrum.size();
//		for (Match match: matchesForOneSpectrum) {
//			score_total += match.getScore();
//		}
//		score_mean = score_total / peptideCount;
//		double difference;
//		for (int i = 0; i < numberOfHistogramBars; i++) {
//			difference = (score_mean - xValues[i]);
//			score_variance +=  difference * difference * histogram[i];
//		}
//		score_standard_deviation = Math.sqrt(score_variance);
//		for (Match match: topMatches) {
//			match.setEValue(getStandardDeviaitonAmount(match.getScore()));
//		}
//	}
//	
//	public double getStandardDeviaitonAmount(double score) {
//		//inverted so lower numbers are better
//		return score_standard_deviation / (score - score_mean);
//	}
//
//}
