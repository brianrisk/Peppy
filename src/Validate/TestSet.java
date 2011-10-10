package Validate;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.Polygon;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

import Peppy.Match;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.ScoringThreadServer;
import Peppy.Spectrum;
import Utilities.U;

public class TestSet {
	
	private String testName;
	private Color color;
	private ArrayList<Spectrum> spectra;
	private ArrayList<MatchContainer> topForwardsTestedMatches = null;
	private ArrayList<Match> topForwardsMatches = null;
	private ArrayList<Match> positiveMatches = new ArrayList<Match>();
	private ArrayList<Match> topReverseMatches = null;
	private ArrayList<MatchContainer> testedMatches = null;
	private int setSize = -1;
	private double averageNumberOfPeaksPerSpectrum;
	
	//statistics
	private int topRankTrueTally = 0;
	private int topRankFalseTally = 0;
	private int totalTrueTally = 0;
	private int trueTallyAtFivePercentError = -1;
	private double percentAtFivePercentError = -1;
	private double eValueAtFivePercentError = -1;
	
	//PR curve
	private double areaUnderPRCurve = 0;
	
	long timeElapsed = 0;
	double milisecondsPerSpectrum;
	
	public TestSet(String testDirectoryName, String testName, ArrayList<Match> matches, Color color) {
		this(testDirectoryName, testName, color);
		positiveMatches = matches;
		cleanMatches();
	}
	

	public TestSet(String testDirectoryName, String testName, Color color) {
		this.testName = testName;
		this.color = color;
		
		//load spectra for this test
		spectra = Spectrum.loadSpectraFromFolder(testDirectoryName + "/" + testName + "/spectra");
		
		setSize = spectra.size();	
		
		averageNumberOfPeaksPerSpectrum = 0;
		for (Spectrum spectrum: spectra) {
			averageNumberOfPeaksPerSpectrum += spectrum.getPeakCount();
		}
		averageNumberOfPeaksPerSpectrum /= spectra.size();
	}
	
	public double getAverageNumberOfPeaksPerSpectrum() {
		return averageNumberOfPeaksPerSpectrum;
	}

	public void findPositiveMatches(ArrayList<Peptide> peptides) {
		/* get the matches and how long it takes */
		long startTimeMilliseconds = System.currentTimeMillis();
//		ArrayList<Match> matches = (new ScoringThreadServer(peptides, spectra)).getMatches();
		ArrayList<Match> matches = Peppy.Peppy.getMatchesWithPeptides(peptides, spectra);
		long stopTimeMilliseconds = System.currentTimeMillis();
		timeElapsed += stopTimeMilliseconds - startTimeMilliseconds;
		
		/* add the matches to the full list */
		if (matches != null) {
			positiveMatches.addAll(matches);
		}
	}

	/**
	 * Here the peptide list is probably from a reverse database.
	 * We are wanting to see what matches we find to a database that does not
	 * have the correct answer.
	 * @param peptides
	 */
	public void findFalsePositiveMatches(ArrayList<Peptide> peptides) {
		Properties.maximumNumberOfMatchesForASpectrum = 1;
		
		//clear E values
		for (Spectrum spectrum: spectra) {
			spectrum.clearEValues();
		}
		//get the matches
		topReverseMatches = (new ScoringThreadServer(peptides, spectra)).getMatches();
		
		//Sort matches by e value	
		Match.setSortParameter(Match.SORT_BY_P_VALUE);
		Collections.sort(topReverseMatches);
		
		//see if any are actually true!
		int trueTally = 0;
		for (Match match: topReverseMatches) {
			MatchContainer mc = new MatchContainer(match);
			if (mc.isTrue()) {
				trueTally++;
				U.p(match);
			}
		}
		if (trueTally > 0) U.p("Some are correct in reverse database! This many: " + trueTally);
		
	}
	

	public void cleanMatches() {
		/* Do a little house cleaning */
//		Peppy.Peppy.assignRankToMatches(positiveMatches);
//		Peppy.Peppy.assignConfidenceValuesToMatches(positiveMatches, spectra);
//		Peppy.Peppy.removeDuplicateMatches(positiveMatches);
		
		/* remove duplicate matches */
//		Match.setSortParameter(Match.SORT_BY_SPECTRUM_ID_THEN_PEPTIDE);
//		Collections.sort(positiveMatches);
//		Match theMatch;
//		Match previousMatch = positiveMatches.get(0);
//		for (int i = 1; i < positiveMatches.size(); i++) {
//			theMatch = positiveMatches.get(i);
//			
//			if (theMatch.getPeptide().equals(previousMatch.getPeptide()) ) {
//				if (theMatch.getSpectrum().getId() == previousMatch.getSpectrum().getId()) {
//					positiveMatches.remove(i);
//					i--;
//				}
//			} else {
//				previousMatch = theMatch;
//			}
//		}
		positiveMatches = Peppy.Peppy.reduceMatchesToOnePerSpectrum(positiveMatches);
		
		/*removing matches with low rank */
//		Match match;
//		int size = positiveMatches.size();
//		for (int i = 0; i < size; i++) {
//			
//			
//			/* get our match */
//			match = positiveMatches.get(i);
//			
//			/* remove matches with low rank */
//			if (match.rank > Properties.maximumNumberOfMatchesForASpectrum) {
//				positiveMatches.remove(i);
//				size--;
//				i--;
//			}
//			
//		}
	}
	
	
	public ArrayList<Match> getPositiveMatches() {
		return positiveMatches;
	}


	public void resetTest() {
		positiveMatches = new ArrayList<Match>();
		topRankTrueTally = 0;
		topRankFalseTally = 0;
		totalTrueTally = 0;
		trueTallyAtFivePercentError = 0;
		percentAtFivePercentError = 0;
		eValueAtFivePercentError = 0;
		areaUnderPRCurve = 0;
	}
	
	/**
	 * This should be called only after our set of matches has been found
	 */
	public void calculateStastics() {		
		//track r time
		milisecondsPerSpectrum = (double) timeElapsed / setSize;

		//sorting by confidence
//		Match.setSortParameter(Match.SORT_BY_RANK_THEN_E_VALUE);
		Match.setSortParameter(Match.SORT_BY_P_VALUE);
//		Match.setSortParameter(Match.SORT_BY_IMP_VALUE);
//		Match.setSortParameter(Match.SORT_BY_SCORE_RATIO);
//		Match.setSortParameter(Match.SORT_BY_P_VALUE);
//		Match.setSortParameter(Match.SORT_BY_RANK_THEN_SCORE);
//		Match.setSortParameter(Match.SORT_BY_SCORE);
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
			if (match.rank == 1) {
				topForwardsMatches.add(match);
				topForwardsTestedMatches.add(new MatchContainer(match));
			}
		}
		//should already be ordered, but what the hell, right?
		Collections.sort(topForwardsMatches);
		Collections.sort(topForwardsTestedMatches);
			
		double highestIMP = 0;
		for (MatchContainer match: topForwardsTestedMatches) {
			if (match.isTrue()) {
				topRankTrueTally++;
				if (match.getMatch().getIMP() > highestIMP) highestIMP = match.getMatch().getIMP();
			} else {
				topRankFalseTally++;
			}
			if ((double) topRankFalseTally / (topRankTrueTally + topRankFalseTally) <= 0.05) {
				trueTallyAtFivePercentError =  topRankTrueTally;
				percentAtFivePercentError = (double) trueTallyAtFivePercentError / setSize;
				eValueAtFivePercentError = match.getEValue();
			}
		}
//		U.p("highestIMP for true found in " + testName + ": " + highestIMP);
		//count total true
		for (MatchContainer match: testedMatches) {
			if (match.isTrue()) totalTrueTally++;
		}
		
		generatePrecisionRecallCurve();
	}


	public BufferedImage generatePrecisionRecallCurve() {
			int width = 500;
			int height = 500;
			//setting up Graphics context
			BufferedImage bufferedImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
			Graphics2D g = bufferedImage.createGraphics();
			g.addRenderingHints(new RenderingHints(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY));
			g.addRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
			g.setColor(Color.white);
			g.fillRect(0,0,width,height);
			
			
			//setting the line color
			g.setColor(color);
			
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
				recall = (double) trueCount / getSetSize();	
//				recall = (double) trueCount / totalTrueTally;	
				
				area += (recall - recallPrevious) * precision;
				
				
				x2 = (int) (recall * width);
				y2 = (int) ((1.0 - precision) * height);
				
				//in case we are moving so little we are not advancing
				if (x1 + 1 >= x2) {
					continue;
				} else {
					x1 = x2_old;
					
	//				U.p(trueCount + ", " + precision);
					g.setColor(color);
	//				g.setStroke(new BasicStroke(2.0f));
					g.drawLine(x1, y1, x2, y2);
					//let's fill in the area under the line, yes?
					Polygon polygon = new Polygon();
					polygon.addPoint(x1, y1);
					polygon.addPoint(x2, y2);
					polygon.addPoint(x2, height);
					polygon.addPoint(x1, height);
					g.setColor(new Color(color.getRed(), color.getGreen(), color.getBlue() ,128));
					g.fillPolygon(polygon);
					
					//updating variables
					x2_old = x2;
					y1 = y2;
				}
				
			}
			//draw final line straight down (just for looks)
			g.setColor(color);
	//		g.setStroke(new BasicStroke(2.0f));
			g.drawLine(x2, y1, x2, height);
			
			//create axis
			int indent = 40;
			int axisSpace = 5;
			BufferedImage axisImage = new BufferedImage(width + indent, height + indent, BufferedImage.TYPE_INT_RGB);
			Graphics2D axisGraphics = axisImage.createGraphics();
			axisGraphics.setFont(new Font("Monospaced", Font.PLAIN, 36));
			FontMetrics fm = axisGraphics.getFontMetrics();
			axisGraphics.setColor(Color.white);
			axisGraphics.fillRect(0, 0, width + indent, height + indent);
			axisGraphics.drawImage(bufferedImage, indent, 0, width, height, null);
			
			//draw axis
			axisGraphics.setColor(Color.black);
			axisGraphics.setStroke(new BasicStroke(3.0f));
			axisGraphics.drawLine(indent, 0, indent, height);
			axisGraphics.drawLine(indent, height, width + indent, height);
			
			//draw axis caps
			int capRadius = 15;
			axisGraphics.drawLine(indent - capRadius, 0, indent + capRadius, 0);
			axisGraphics.drawLine(indent + width, height - capRadius, indent + width, height + capRadius);
			
			//draw labels
			axisGraphics.drawString("0", indent - fm.stringWidth("0") - axisSpace, height + fm.getHeight());
			axisGraphics.drawString("1", indent - fm.stringWidth("1") - axisSpace, fm.getHeight() + axisSpace);
			axisGraphics.drawString("1", width + indent - fm.stringWidth("1") - axisSpace, height + fm.getHeight());
			
			areaUnderPRCurve = area;
			return axisImage;
			
		}

	public ArrayList<Peptide> loadCorrectPeptides() {
		ArrayList<Peptide> correctPeptides = new ArrayList<Peptide>();
		for(Spectrum spectrum: spectra) {
			//find the file for the correct peptide
			File spectrumFile = spectrum.getFile();
			File testFolder = spectrumFile.getParentFile().getParentFile();
			File peptideFolder = new File(testFolder, "peptides");
			File peptideFile = new File(peptideFolder, spectrumFile.getName());
			
			//load in the correct peptide string
			boolean validPeptideFile = true;
			String correctAcidSequence = "";
			try {
				BufferedReader br = new BufferedReader(new FileReader(peptideFile));
				//read the first line;
				correctAcidSequence = br.readLine();
				//close;
				br.close();
			} catch (FileNotFoundException e) {
				validPeptideFile = false;
				e.printStackTrace();
			} catch (IOException e) {
				validPeptideFile = false;
				e.printStackTrace();
			}
			
			//testing that we've got a valid peptide file
			if (correctAcidSequence == null) {
				validPeptideFile = false;
			}
			correctAcidSequence = correctAcidSequence.trim();
			if (correctAcidSequence.equals("")) {
				validPeptideFile = false;
			}
			//adding to the array list
			if (validPeptideFile) {
				correctPeptides.add(new Peptide(correctAcidSequence));
			}
		}
		return correctPeptides;
	}
	
	@SuppressWarnings("unused")
	private ArrayList<Match> loadCorrectMatches() {
		ArrayList<Match> correctMatches = new ArrayList<Match>();
		for(Spectrum spectrum: spectra) {
			//find the file for the correct peptide
			File spectrumFile = spectrum.getFile();
			File testFolder = spectrumFile.getParentFile().getParentFile();
			File peptideFolder = new File(testFolder, "peptides");
			File peptideFile = new File(peptideFolder, spectrumFile.getName());
			
			//load in the correct peptide string
			boolean validPeptideFile = true;
			String correctAcidSequence = "";
			try {
				BufferedReader br = new BufferedReader(new FileReader(peptideFile));
				//read the first line;
				correctAcidSequence = br.readLine();
				//close;
				br.close();
			} catch (FileNotFoundException e) {
				validPeptideFile = false;
				e.printStackTrace();
			} catch (IOException e) {
				validPeptideFile = false;
				e.printStackTrace();
			}
			
			//testing that we've got a valid peptide file
			if (correctAcidSequence == null) {
				validPeptideFile = false;
			}
			correctAcidSequence = correctAcidSequence.trim().toUpperCase();
			if (correctAcidSequence.equals("")) {
				validPeptideFile = false;
			}
			
			//testing we've got a valid peptide
			boolean validPeptide = true;
			Peptide peptide = new Peptide(correctAcidSequence);
			if (peptide.getMass() < 0) {
				validPeptide = false;
				U.p(correctAcidSequence);
				U.p(testName);
				U.p(spectrumFile.getName());
			}
			//adding to the array list
			if (validPeptideFile && validPeptide) {
				correctMatches.add(Properties.matchConstructor.createMatch(spectrum, peptide));
			}
			
//			if (spectrum.getFile().getName().equals("T10707_Well_H13_1768.77_19185.mgf..pkl")) {
//				U.p ("Valid? " + validPeptideFile);
//			}
		}
		return correctMatches;
	}
	
	
	public String getName() {
		return testName;
	}

	public int getTrueTally() {
		return topRankTrueTally;
	}

	public int getFalseTally() {
		return topRankFalseTally;
	}

	public int getTrueTallyAtFivePercentError() {
		return trueTallyAtFivePercentError;
	}

	public double getEValueAtFivePercentError() {
		return eValueAtFivePercentError;
	}

	public double getPercentAtFivePercentError() {
		return percentAtFivePercentError;
	}

	public double getMilisecondsPerSpectrum() {
		return milisecondsPerSpectrum;
	}

	public double getAreaUnderPRCurve() {
		return areaUnderPRCurve;
	}

	public long getTimeElapsed() {
		return timeElapsed;
	}

	public int getSetSize() {
		return setSize;
	}

	public ArrayList<MatchContainer> getTestedMatches() {
		return testedMatches;
	}
	
	public ArrayList<MatchContainer> getTopForwardsTestedMatches() {
		return topForwardsTestedMatches;
	}

	public double getEValueAtPercentForwards(double percent) {
		Match.setSortParameter(Match.SORT_BY_P_VALUE);
		Collections.sort(topForwardsMatches);
		int level = (int) (topForwardsMatches.size() * percent);
		return topForwardsMatches.get(level).getEValue();
	}
	
	public double getEValueAtPercentReverse(double percent) {
		Match.setSortParameter(Match.SORT_BY_P_VALUE);
		Collections.sort(topReverseMatches);
		int level = (int) (topReverseMatches.size() * percent);
		return topReverseMatches.get(level).getEValue();
	}


}
