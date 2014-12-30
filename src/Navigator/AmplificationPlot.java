package Navigator;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;

import javax.imageio.ImageIO;

import Experimental.AmplifiedPeptide;
import Peppy.U;

/**
 * I am interested in how amplification as measured by aCGH affects expression levels.
 * 
 * The idea is to search the W2 and W16 data on the personal proteome files that contain
 * the aCGH amplification scores.  We will then create hashtables of peptides,integers where the
 * integer is the count for how many peptides were found.  We will do this for both W2 and W16.
 * 
 * The hashtable will only contain peptides that are in some area of amplification
 * 
 * 
 * @author Brian Risk
 *
 */
public class AmplificationPlot {
	
	/* how much do we want to split these bins?  into quarters? units of 1?  This is the divisor */
	static int divisor = 5;
	static int maxValue = 10;
	
	
	static int whim2Total = 0;
	static int whim16Total = 0;
	static Hashtable<String, AmplifiedPeptide> whim2Hash = new Hashtable<String, AmplifiedPeptide>();
	static Hashtable<String, AmplifiedPeptide> whim16Hash = new Hashtable<String, AmplifiedPeptide>();
	static ArrayList<Point2D.Double> peptidePoints = new ArrayList<Point2D.Double>();
	static Hashtable<Double, ArrayList<Double>> bins = new Hashtable<Double, ArrayList<Double>>();
	
	public static void main(String [] args) {
		U.p("finding amplifciation plot");
		U.p("loading files...");
		
		ArrayList<File> w2Files = new ArrayList<File>();
//		w2Files.add(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/Amp-UNC-WHIM2-Ellis043/1 WHIM2 - refPersonalProteomeOutput.fasta/report.txt"));
		w2Files.add(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/Amp-WHIM2-Ellis043/1 WHIM2 - refPersonalProteomeOutput.fasta/report.txt"));
//		w2Files.add(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/Amp-WHIM2-Ellis041/1 WHIM2 - refPersonalProteomeOutput.fasta/report.txt"));
//		w2Files.add(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/Amp-WHIM2-Ellis033/1 WHIM2 - refPersonalProteomeOutput.fasta/report.txt"));
		
		ArrayList<File> w16Files = new ArrayList<File>();
//		w16Files.add(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/Amp-UNC-WHIM16-Ellis043/1 WHIM16 - refPersonalProteomeOutput.fasta/report.txt"));
		w16Files.add(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/Amp-WHIM16-Ellis043/1 WHIM16 - refPersonalProteomeOutput.fasta/report.txt"));
//		w16Files.add(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/Amp-WHIM16-Ellis041/1 WHIM16 - refPersonalProteomeOutput.fasta/report.txt"));
//		w16Files.add(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/Amp-WHIM16-Ellis033/1 WHIM16 - refPersonalProteomeOutput.fasta/report.txt"));
		
		
		
		
		
		for (File file: w16Files) {
			whim16Total += addMatches(file, whim16Hash);
		}
		for (File file: w2Files) {
			whim2Total += addMatches(file, whim2Hash);
		}
		
		/* which hash comes first is arbitrary */
		createPeptidePoints(whim2Hash, whim16Hash);
		createBins();
		createReport();
		
		drawGrid(100,100);
		
		U.p("done");
	}
	
	private static int addMatches(File matchesFile, Hashtable<String, AmplifiedPeptide> hashtable) {
		ArrayList<Match> fileMatches = Match.loadMatches(matchesFile);
		for (Match match: fileMatches) {
			
			AmplifiedPeptide amp = new AmplifiedPeptide(match);
			if (!amp.isValid()) continue;
			AmplifiedPeptide storedAmp = hashtable.get(amp.getPeptide());
			if (storedAmp == null) {
				hashtable.put(amp.getPeptide(), amp);
			} else {
				storedAmp.increment();
			}
		}
		return fileMatches.size();
	}
	
	/**
	 * Finds all unique peptides where in one sample there is amplification, and in the other there isn't.
	 * creates a data point for each of these peptides.
	 * @param hashA
	 * @param hashB
	 */
	private static void createPeptidePoints(Hashtable<String, AmplifiedPeptide> hashA, Hashtable<String, AmplifiedPeptide> hashB) {
		ArrayList<AmplifiedPeptide> peptidesA = new ArrayList<AmplifiedPeptide>(hashA.values());
		for (AmplifiedPeptide ampA: peptidesA) {
			AmplifiedPeptide ampB = hashB.get(ampA.getPeptide());
			
			/* we only want cases where both samples have that peptide and only one falls in an amplified region */
			if (ampB == null) continue;
			boolean aIsAmplified = (ampA.getScore() == 0);
			boolean bIsAmplified = (ampB.getScore() == 0);
			if (!aIsAmplified && !bIsAmplified) continue;
			if (aIsAmplified && bIsAmplified) continue;
			if (ampA.getCount() < 10 || ampB.getCount() < 10) continue;
			
			/* x is the amp score */
			double xValue = 0;
			
			/* y is the proportional difference */
			double yValue = 0;
			
			/*NOTE this violates the arbitrary order */
			double normalizedA = (double) ampA.getCount() / (double) whim2Total;
			double normalizedB = (double) ampB.getCount() / (double) whim16Total;
			
			/* if ampA's score is zero, that means ampB is the one with the amplification going on */
			if (aIsAmplified) {
				xValue = ampB.getScore();
				if (ampA.getCount() > ampB.getCount()) {
					yValue = (double) normalizedA / (double) normalizedB * -1.0 + 1.0;
				} else {
					yValue = (double) normalizedB / (double) normalizedA - 1.0;
				}
			} else {
				xValue = ampA.getScore();
				if (ampA.getCount() > ampB.getCount()) {
					yValue = (double) normalizedA / (double) normalizedB - 1.0;
				} else {
					yValue = (double) normalizedB / (double) normalizedA * -1.0 + 1.0;
				}
			}	
			
			if(xValue < -0.5 && yValue > 5) U.p(ampA.getPeptide() + ", x: " + xValue + ", y: " + yValue + ", W2: " + aIsAmplified);
			
			peptidePoints.add(new Point2D.Double(xValue, yValue));
		}
	}
	
	/**
	 * we want to see the distribution of things within a range of amplification.
	 * This bins up the data.
	 */
	public static void createBins() {
		for (Point2D.Double point: peptidePoints) {
			double bin = getBin(point.x);
			ArrayList<Double> binEntries = bins.get(bin);
			if (binEntries == null) {
				binEntries = new ArrayList<Double>();
				bins.put(bin, binEntries);
			}
			binEntries.add(point.y);
		}
	}
	
		
	private static void createReport() {
		
		/* amplification bin plot.txt */
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("amplification bin plot.txt")));
			ArrayList<Double> binKeys = new ArrayList<Double>(bins.keySet());
			for (Double binKey: binKeys) {
				ArrayList<Double> binValues = bins.get(binKey);
				Point2D.Double meanSD = getMeanAndStandardDeviation(binValues);
				pw.println(binKey + "\t" + meanSD.x + "\t" + meanSD.y + "\t" + binValues.size());
			}
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		/* amplification point plot.txt */
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("amplification point plot.txt")));
			for (Point2D.Double point: peptidePoints) {
				pw.println(point.x + "\t" + point.y);
			}
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		/* create the heat map */
		int matrixWidth = maxValue * divisor * 2;
		U.p("matrixWidth: " + matrixWidth);
		int[][] countArray = new int[matrixWidth][matrixWidth];
		int xIndex, yIndex;
		for (Point2D.Double point: peptidePoints) {
			xIndex = (int) (divisor * (point.x +  maxValue));
			yIndex = (int) (divisor * (point.y +  maxValue));
			if (xIndex >= matrixWidth) continue;
			if (xIndex < 0) continue;
			if (yIndex >= matrixWidth) continue;
			if (yIndex < 0) continue;
			countArray[xIndex][yIndex]++;
		}
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("amplificaiton grid values.txt")));
			for (yIndex = 0; yIndex < matrixWidth; yIndex++) {
				StringBuffer line = new StringBuffer();
				for(xIndex = 0; xIndex < matrixWidth; xIndex++) {
					line.append(countArray[xIndex][yIndex] + "\t");
				}
				pw.println(line);
			}
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	private static void drawGrid(int rows, int cols) {
		/* load the grid */
		double [][] results = new double [rows][cols];
		
		/* load in the results matrix */
		try {
			BufferedReader br = new BufferedReader(new FileReader("amplificaiton grid values.txt"));
			String line = br.readLine();
			int rowIndex = 0;
			while (line != null) {
				String [] chunks = line.split("\t");
				for (int colIndex = 0; colIndex < chunks.length; colIndex++) {
					results[rowIndex][colIndex] = Integer.parseInt(chunks[colIndex]);
				}
				
				line = br.readLine();
				rowIndex++;
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		/* take the log! */
		double maxLog = 0;
		for (int rowIndex = 0; rowIndex < rows; rowIndex++) {
			for (int colIndex = 0; colIndex < cols; colIndex++) {
				if (results[rowIndex][colIndex] != 0) {
					results[rowIndex][colIndex] = Math.log(results[rowIndex][colIndex]);
					results[rowIndex][colIndex] += 1;
					if (results[rowIndex][colIndex] > maxLog) maxLog = results[rowIndex][colIndex];
				}
			}
		}
		
		/* create percent matrix */
		double [][] percentages = new double[rows][cols];
		
		
		for (int rowIndex = 0; rowIndex < rows; rowIndex++) {
			double rowSum = 0;
			double rowMax = 0;
			for (int colIndex = 0; colIndex < cols; colIndex++) {
				rowSum += results[rowIndex][colIndex];
				if (results[rowIndex][colIndex] > rowMax) rowMax = results[rowIndex][colIndex];
			}
			double normalizingFactor = (double) rowSum / rowMax;
			for (int colIndex = 0; colIndex < cols; colIndex++) {
				double percentage = 0;
				if (rowSum != 0) {
					percentage = (double) results[rowIndex][colIndex] / maxLog;
//					percentage *= normalizingFactor;
					if (percentage > 1) percentage = 1;
				}
				percentages[rowIndex][colIndex] = percentage;
			}
		}
		
		int cellWidth = 10;
		int cellHeight = cellWidth;
		int verticalScale = 1;
		int imageWidth = cellWidth * cols;
		int imageHeight = cellHeight * rows / verticalScale;
		BufferedImage bdest = new BufferedImage(imageWidth, imageHeight, BufferedImage.TYPE_INT_RGB);
		Graphics2D g = bdest.createGraphics();
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_SPEED));
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF));
		for (int rowIndex = 0; rowIndex < rows; rowIndex += verticalScale) {
			for (int colIndex = 0; colIndex < cols; colIndex++) {
				Color color;
				if (percentages[rowIndex][colIndex] == 0) {
					color = Color.black;
				} else {
					color = Color.getHSBColor((float) percentages[rowIndex][colIndex], 1, 1);
					int colorValue = (int) (percentages[rowIndex][colIndex] * 255.0);
					color = new Color(colorValue, colorValue, colorValue);
				}
				int x = colIndex * cellWidth;
				int y = imageHeight - (rowIndex * cellHeight / verticalScale);
				g.setColor(color);
				g.fillRect(x, y, cellWidth, cellHeight);
			}
		}
		/* draw axis */
		g.setColor(Color.magenta);
		int halfWidth = (cellWidth * cols) / 2;
		g.drawLine(0, halfWidth, imageWidth, halfWidth);
		g.drawLine(halfWidth, 0, halfWidth, imageHeight);
		/* tics */
		int ticTop = halfWidth - cellWidth;
		int ticBottom = halfWidth + cellWidth;
		int ticLoc;
		for (int i = 0; i < maxValue * 2; i++) {
			ticLoc = i * divisor * cellWidth;
			g.drawLine(ticLoc, ticTop, ticLoc, ticBottom);
			g.drawString("" + (i - maxValue), ticLoc, ticTop);
			g.drawLine(ticTop, ticLoc, ticBottom, ticLoc);
			g.drawString("" + (-i + maxValue), ticBottom, ticLoc);
		}
		try {
			ImageIO.write(bdest,"PNG",new File("probailityGrid.png"));
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private static double getBin(double value){
		return  Math.ceil(value * (double) divisor) / (double) divisor;
	}
	
	
	private static Point2D.Double getMeanAndStandardDeviation(ArrayList<Double> values) {
		double mean = 0;
		for (Double value: values) {
			mean += value;
		}
		mean = mean / values.size();
		
		double standardDeviation = 0;
		double difference;
		for (Double value: values) {
			difference = mean - value;
			standardDeviation += difference * difference;
		}
		standardDeviation /= values.size();
		standardDeviation = Math.sqrt(standardDeviation);
		
		return new Point2D.Double(mean, standardDeviation);		
	}
	

}
