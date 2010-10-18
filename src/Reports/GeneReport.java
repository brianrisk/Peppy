package Reports;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import Peppy.Match;
import Peppy.Properties;
import Peppy.Sequence;
import Peppy.SequenceRegion;
import Utilities.U;

public class GeneReport {
	Sequence sequence;
	ArrayList<Match> matches;
	int geneWindowSize;
	int [] indexScores;
//	int [] scoresForWindow;
	int [] windowScores;
	
	//E-Values:  allocating histogram variables
	private final int numberOfHistogramBars = 100;
	private int [] eValueHistogram = new int[numberOfHistogramBars];
	private double [] scoreProbabilities = new double[numberOfHistogramBars];
	private double [] survivability = new double[numberOfHistogramBars];
	private double [] xValues = new double[numberOfHistogramBars];
	private double barWidth;
	@SuppressWarnings("unused")
	private double m, b;
	
	public GeneReport(Sequence sequence, ArrayList<Match> matches, int geneRegionSize) {
		this.sequence = sequence;
		this.matches = matches;
		this.geneWindowSize = geneRegionSize;
		indexScores = new int[sequence.getSequenceLength()];
		for (int i = 0; i < indexScores.length ; i++) {
			indexScores[i] = 0;
		}
		populateIndexScores();
		calculateWindowScores();
		calculateEValues();
	}
	
	public void populateIndexScores() {
		int score, peptideIndex;
		for (Match match: matches) {
			score = (int) Math.round(Math.abs(Math.log(match.getEValue())));
			peptideIndex = match.getPeptide().getStartIndex();
			indexScores[peptideIndex] += score;	
		}
	}
	
	public void calculateWindowScores() {
		int windowScoresLength = indexScores.length - geneWindowSize;
		windowScores = new int[windowScoresLength];
		for (int i = 0; i < windowScoresLength; i++) {
			windowScores[i] = 0;
		}
		
		//find value for first window
		int windowTally = 0;
		for (int i = 0; i < geneWindowSize; i++) {
			windowTally += indexScores[i];
		}
		windowScores[0] = windowTally;
		
		//find the rest of the windows
		for (int i = 1; i < windowScoresLength; i++) {
			windowTally -= indexScores[i - 1];
			windowTally += indexScores[i + geneWindowSize - 1];
			windowScores[i] = windowTally;
		}
	}
	
	public void insertSequenceRegionReport(ArrayList<Match> theseMatches, File sequenceReportDirectory, PrintWriter pw) throws IOException {
		//initialize our variables
//		int numberOfTopRegions = 10;
		int histogramHeight = 200;
//		int sequenceLength = sequence.getSequenceLength();
		
		int bestWindowIndex = 0; 
		int bestWindowValue = 0;
		for (int i = 0; i < windowScores.length; i++) {
			if (windowScores[i] > bestWindowValue) {
				bestWindowValue = windowScores[i];
				bestWindowIndex = i;
			}
		}
		int bestWindowIndexTop = bestWindowIndex + geneWindowSize - 1;
		
		
		//fill our regions with matches
		int matchIndex;
		SequenceRegion topSequenceRegion = new SequenceRegion(bestWindowIndex, geneWindowSize);
		for (Match match: theseMatches) {
			matchIndex = match.getPeptide().getStartIndex();
			if (matchIndex >= bestWindowIndex && matchIndex <= bestWindowIndexTop) {
				topSequenceRegion.addHit(match);
			}
		}

		
		//create visualizations and report of the top regions
		File regionDirectory = new File(sequenceReportDirectory, "regions" + sequence.getId());
		regionDirectory.mkdirs();
		pw.println("<table border=1>");
//		for (int i = 0; i < numberOfTopRegions; i++) {
//			region = regions.get(i);
			int i = 0;
			pw.println();
			pw.print("<tr>");
			
			pw.print("<td>");
			pw.print("<a href=\"" + regionDirectory.getName() + "/" + i + ".html\">");
			pw.print(topSequenceRegion.getStartIndex());
			pw.print("</a>");
			pw.print("</td>");
			
			pw.print("<td>");
			pw.print(topSequenceRegion.getScore());
			pw.print("</td>");
			
			pw.print("<td>");
			pw.print(topSequenceRegion.getMaxBar());
			pw.print("</td>");
			
			pw.print("<td>");
			pw.print("<img src=\"" + regionDirectory.getName() + "/" + i + ".jpg\">");
			pw.print("</td>");
			
			pw.print("</tr>");
			
			//create the images and html for this region
			HistogramVisualizer.drawShadedHistogram(topSequenceRegion.getHistogram(), histogramHeight, topSequenceRegion.getMaxBar(), new File(regionDirectory, i + ".jpg"));
//		}
		pw.println("</table>");
	}
	
	/**
	 * find expected value (a.k.a. "e value") for top matches
	 * 
	 * This method assumes that matchesForOneSpectrum is already sorted from highest score to lowest.
	 * Calculates e values for each of the top matches
	 * @param matchesForOneSpectrum
	 * @param topMatches
	 */
	public void calculateEValues() {
		int highScore = 0;
		int lowScore = Integer.MAX_VALUE;
		for (int i = 0; i < windowScores.length; i++) {
			if (windowScores[i] < 20) continue;
			if (windowScores[i] > highScore) highScore = windowScores[i]; 
			if (windowScores[i] < lowScore) lowScore = windowScores[i]; 
		}

		barWidth = (highScore - lowScore) / numberOfHistogramBars;
		
		//initializing histograms and xValues
		for (int i = 0; i < numberOfHistogramBars; i++) {
			eValueHistogram[i] = 0;
			xValues[i] = lowScore + (i * barWidth);
		}


		//populate the histogram
		int bin;
		for (int i = 0; i < windowScores.length; i++) {
			if (windowScores[i] < 20) continue;
			bin = (int) Math.floor((windowScores[i] - lowScore) / barWidth);
			if (bin < numberOfHistogramBars) {
				eValueHistogram[bin]++;
			} else {
				eValueHistogram[numberOfHistogramBars - 1]++;
			}
//			scoreSummation += windowScores[i];
		}
		

		File histogramFile = new File(Properties.reportDirectory, sequence.getId() + "-hist.jpg");
		try {
			HistogramVisualizer.drawHistogram(eValueHistogram, 300, 300, histogramFile);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		//find score probabilities
		for (int i = 0; i < numberOfHistogramBars; i++) {
			scoreProbabilities[i] = (double) eValueHistogram[i] / windowScores.length;
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
			if (indexScores[topIndex] == 0) break;
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
		m = U.calculateM(xValues, survivability, chopIndex, topIndex);
		b = U.calculateB(xValues, survivability, chopIndex, topIndex, m);
		
		//using our m and b to derive e values for all top matches
//		double eValue;
//		for (Match match: topMatches) {
//			eValue = m * match.getScore() + b;
//			eValue = Math.exp(eValue);
//			eValue *= windowScores.length;
//			//setting to -1 if eValue is Nan
//			if (eValue <= 1 == eValue >= 1) eValue = Double.MAX_VALUE;
//			match.setEValue(eValue);
//		}
	}
	
	public int [] getEValueHistogram() {return eValueHistogram;}

}
