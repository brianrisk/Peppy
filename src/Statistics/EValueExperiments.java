package Statistics;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Polygon;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;

import javax.imageio.ImageIO;

import Peppy.Match;
import Peppy.Peppy;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.ScoringThreadServer;
import Peppy.Sequence;
import Peppy.Spectrum;
import Utilities.U;
import Validate.MatchContainer;
import Validate.ReliabilityTester;

public class EValueExperiments {
	
	private final int numberOfHistogramBars = 100;
	
	private static String testName = "human";
	private ArrayList<Spectrum> spectra;
	private ArrayList<MatchContainer> topForwardsTestedMatches = null;
	private ArrayList<Match> topForwardsMatches = null;
	private ArrayList<Match> positiveMatches = null;
	private ArrayList<MatchContainer> testedMatches = null;
	private int setSize = -1;
	
	//statistics
	private int trueTally = 0;
	private int falseTally = 0;
	private int trueTallyAtOnePercentError = -1;
	private double percentAtOnePercentError = -1;
	private double eValueAtOnePercentError = -1;
	
	//PR curve
	private double areaUnderPRCurve = 0;
	
	public static void main(String args[]) {
		setProperties();
		saveHistograms();
//		new EValueExperiments();
	}
	
	public static void setProperties() {
		Properties.spectrumToPeptideMassError = 2;
		Properties.isSequenceFileDNA = false;
		Properties.spectraDirectoryOrFile = new File("/Users/risk2/PeppyOverflow/tests/" + testName + "/spectra");
		Properties.sequenceDirectoryOrFile = new File("/Users/risk2/PeppyOverflow/tests/databases/uniprot_sprot.fasta");
		//we'd prefer not to have duplicate matches -- especially for the correct ones
		Properties.reduceDuplicateMatches = true;
	}
	
	public static void saveHistograms() {
		//Load our spectra
		U.p("loading spectra...");
		ArrayList<Spectrum> spectra = Spectrum.loadSpectra();
		Collections.sort(spectra);
		U.p("loaded " +spectra.size() + " spectra.");
		
		//Get references to our sequence files -- no nucleotide data is loaded at this point
		ArrayList<Sequence> sequences = Sequence.loadSequences(Properties.sequenceDirectoryOrFile);
		
		//initialize our ArrayList of matches
		ArrayList<Match> matches = new ArrayList<Match>();
		
		//loop through each sequence in the sequences ArrayList
		for (Sequence sequence: sequences) {		
			matches.addAll(Peppy.getMatches(sequence, spectra));
		}
		
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("evalue-tests/histograms.txt")));
			for (Spectrum spectrum: spectra) {
				int [] histogram = spectrum.getEValueCalculator().getHistogram();
				StringBuffer sb = new StringBuffer();
				for (int i = 0; i < histogram.length; i++) {
					sb.append(histogram[i]);
					sb.append('\t');
				}
				pw.println(sb.toString());
			}
			pw.flush();
			pw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		U.p("done");
	}

	public EValueExperiments() {
		ArrayList<Peptide> peptides = ReliabilityTester.loadHighScoringPeptides(testName);
		findPositiveMatches(peptides);
		printStatistics();
		drawHistograms();
	}
	
	private void drawHistograms() {
		
	}
	
	public static void drawHistogram(int [] histogram, int height, int width, File dest) throws IOException {
		//setting up Graphics context
		BufferedImage bdest = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		Graphics2D g = bdest.createGraphics();
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_SPEED));
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF));
		g.setColor(Color.white);
		g.fillRect(0,0,width,height);
		
		//find the maximum value
		int max = 0;
		for (int i = 0; i < histogram.length; i++) {
			if (histogram[i] > max) max = histogram[i];
		}
		
		//draw the bars
		g.setColor(Color.darkGray);
		g.drawRect(0,0,width - 1, height - 1);
		double scaleFactor = (double) height / max;
		for (int i = 0; i < histogram.length; i++) {
			int x1 = (int) ((double) i * width / histogram.length);
			int barWidth = (int) ((double) (i + 1.0) * width / histogram.length) - x1;
			if (barWidth < 1) barWidth = 1;
			int barHeight = (int) (histogram[i] * scaleFactor);
			//account for tiny bars that might not be displayed
			if (barHeight <= 1 && histogram[i] > 0) {
				barHeight = 5;
				g.setColor(Color.lightGray);
			} else {
				g.setColor(Color.black);
			}
			g.fillRect(x1, height - barHeight, barWidth, barHeight);
		}
		
		//outline
		g.setColor(Color.black);
		g.drawRect(0,0,width - 1, height - 1);
		
		//write
		ImageIO.write(bdest,"JPG",dest);
	}
	
	public void findPositiveMatches(ArrayList<Peptide> peptides) {

		//load spectra for this test
		spectra = Spectrum.loadSpectra();
		Collections.sort(spectra);
		setSize = spectra.size();
		U.p("number of spectra: " + setSize);
		
		//get the matches
		positiveMatches = (new ScoringThreadServer(peptides, spectra, null)).getMatches();
		
		//load histograms; recalculate EValues
		try {
			BufferedReader br = new BufferedReader(new FileReader("evalue-tests/histograms.txt"));
			String line = br.readLine();
			//there should be a line for each spectrum
			for (Spectrum spectrum: spectra) {
				//get the histogram
				String [] chunks = line.split("\t");
				int [] histogram = new int[numberOfHistogramBars];
				for (int i = 0; i < numberOfHistogramBars; i++) {
					histogram[i] = Integer.valueOf(chunks[i].trim());
				}
				//set the histogram for the spectrum
				spectrum.getEValueCalculator().setHistogram(histogram);
				
				//read the next line
				line = br.readLine();
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		U.p("calculating final e values");
		for (Match match: positiveMatches) {
			U.p(match.calculateEValue());
		}
		

		//Sort matches by e value	
		Match.setSortParameter(Match.SORT_BY_E_VALUE);
		Collections.sort(positiveMatches);
		
		//See which of the positive matches are true
		testedMatches = new ArrayList<MatchContainer>();
		for (Match match: positiveMatches) {
			if (match == null) U.p("null match?  What?");
			testedMatches.add(new MatchContainer(match));
		}
		//Sort MatchContainer by e value (default)	
		Collections.sort(testedMatches);
		
		//find the #1 ranked match for each spectrum
		topForwardsMatches = new ArrayList<Match>();
		topForwardsTestedMatches = new ArrayList<MatchContainer>();
		for (Match match: positiveMatches) {
			if (match.getRank() == 0) {
				topForwardsMatches.add(match);
				topForwardsTestedMatches.add(new MatchContainer(match));
			}
		}
		//should already be ordered, but what the hell, right?
		Collections.sort(topForwardsMatches);
		Collections.sort(topForwardsTestedMatches);
			
		calculateStastics();
	}
	
	/**
	 * This should be called only after our set of matches has been found
	 */
	private void calculateStastics() {
		//stats for tested matches
		boolean onePercenThresholdHasBeenReached = false;
		for (MatchContainer match: topForwardsTestedMatches) {
			if (match.isTrue()) {
				trueTally++;
			} else {
				falseTally++;
			}
			if (!onePercenThresholdHasBeenReached) {
				if ((double) falseTally / setSize >= 0.01) {
					onePercenThresholdHasBeenReached = true;
					trueTallyAtOnePercentError =  trueTally;
					percentAtOnePercentError = (double) trueTallyAtOnePercentError / setSize;
					eValueAtOnePercentError = match.getEValue();
				}
			}
		}
		
		generatePrecisionRecallCurve();
	}
	
	public void printStatistics() {
		U.p("percentAtOnePercentError:" + percentAtOnePercentError);
		U.p("eValueAtOnePercentError:" + eValueAtOnePercentError);
		U.p("areaUnderPRCurve:" + areaUnderPRCurve);
	}
	
	private void generatePrecisionRecallCurve() {
		int width = 500;
		int height = 500;
		//setting up Graphics context
		BufferedImage bdest = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		Graphics2D g = bdest.createGraphics();
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY));
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
		g.setColor(Color.white);
		g.fillRect(0,0,width,height);
		
		//setting the line color
		g.setColor(Color.red);
		
		int x1 = 0;
		int x2 = 0;
		int x2_old = 0;
		int y1 = 0;
		int y2 = 0;
		int trueCount = 0;
		double precision = 0; 
		double recall = 0;
		double recallPrevious = 0;
		double area = 0;
		for(int i = 0; i < testedMatches.size(); i++) {
			MatchContainer match = testedMatches.get(i);
			if (match.isTrue()) {
				trueCount++;
			}
			
			
			precision = (double) trueCount / (i + 1);
			recallPrevious = recall;
			recall = (double) trueCount / setSize;	
			
			area += (recall - recallPrevious) * precision;
			
			
			x2 = (int) (recall * width);
			y2 = (int) ((1.0 - precision) * height);
			
			//in case we are moving so little we are not advancing
			if (x1 + 1 >= x2) {
				continue;
			} else {
				x1 = x2_old;
				
//				U.p(trueCount + ", " + precision);
				g.setColor(Color.red);
//				g.setStroke(new BasicStroke(2.0f));
				g.drawLine(x1, y1, x2, y2);
				//let's fill in the area under the line, yes?
				Polygon polygon = new Polygon();
				polygon.addPoint(x1, y1);
				polygon.addPoint(x2, y2);
				polygon.addPoint(x2, height);
				polygon.addPoint(x1, height);
				g.setColor(new Color(255,0,0,128));
				g.fillPolygon(polygon);
				
				//updating variables
				x2_old = x2;
				y1 = y2;
			}
			
		}
		//draw final line straight down (just for looks)
		g.setColor(Color.red);
//		g.setStroke(new BasicStroke(2.0f));
		g.drawLine(x2, y1, x2, height);
		
		areaUnderPRCurve = area;
		
//		g.drawLine(x2, y2, width, height);
		
		try {
			File testDirectory = new File(Properties.validationDirectory, testName);
			testDirectory.mkdirs();
			File imageFile = new File("evalue-tests/precision-recall.jpg");
			ImageIO.write(bdest,"JPG",imageFile);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
