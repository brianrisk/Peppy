package Validate;

import java.awt.Color;
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

import javax.imageio.ImageIO;

import Peppy.Match;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.ScoringThreadServer;
import Peppy.Spectrum;
import Utilities.U;

public class TestSet {
	
	private String testName;
	private ArrayList<Spectrum> spectra;
	private ArrayList<MatchContainer> topForwardsTestedMatches = null;
	private ArrayList<Match> topForwardsMatches = null;
	private ArrayList<Match> positiveMatches = null;
	private ArrayList<Match> topReverseMatches = null;
	private ArrayList<Match> correctMatches = null;
	private ArrayList<MatchContainer> testedMatches = null;
	private int setSize = -1;
	
	//statistics
	private int topRankTrueTally = 0;
	private int topRankFalseTally = 0;
	private int totalTrueTally = 0;
	private int trueTallyAtOnePercentError = -1;
	private double percentAtOnePercentError = -1;
	private double eValueAtOnePercentError = -1;
	
	//PR curve
	private double areaUnderPRCurve = 0;
	private String fileNameForPRCurve = "";
	
	long timeElapsed;
	String timeToComplete = "";
	double milisecondsPerSpectrum;
	
	public TestSet(String testName) {
		this.testName = testName;
		
		//load spectra for this test
		spectra = Spectrum.loadSpectraFromFolder("/Users/risk2/PeppyOverflow/tests/" + testName + "/spectra");
		setSize = spectra.size();
		
//		//load "correct" matches as defined by the test set
//		correctMatches = loadCorrectMatches();
//		
//		//use only spectra that don't have a big difference between precursor and the predicted protein mass
//		ArrayList<Spectrum> reducedSpectra = new ArrayList<Spectrum>();
//		for (Match match: correctMatches) {
//			double difference = match.getSpectrum().getPrecursorMass() - match.getPeptide().getMass();
//			if (Math.abs(difference) < Properties.peptideMassThreshold) {
//				reducedSpectra.add(match.getSpectrum());
//			}
//		}
//		spectra = reducedSpectra;
//		setSize = spectra.size();
//		
//		//add all the correct peptides to the peptide database
//		peptides.addAll(loadCorrectPeptides());
//		Collections.sort(peptides);		
	}
	
	public void findPositiveMatches(ArrayList<Peptide> peptides) {
//		//first remove spectra which do not represent peptides in our given database
//		correctMatches = loadCorrectMatches();
//		for (Match match: correctMatches) {
//			Peptide peptide = match.getPeptide();
//			if (GenerateValidationReport.isPeptidePresentInList(peptide, peptides) == -1) {
//				for (int i = 0; i < spectra.size(); i++) {
//					if (spectra.get(i).getId() == match.getSpectrum().getId()) {
//						spectra.remove(i);
//						i--;
//					}	
//				}
//			}
//		}
//		setSize = spectra.size();
		
		
		//get the matches
		long startTimeMilliseconds = System.currentTimeMillis();
		positiveMatches = (new ScoringThreadServer(peptides, spectra, null)).getMatches();
		long stopTimeMilliseconds = System.currentTimeMillis();
		timeElapsed = stopTimeMilliseconds - startTimeMilliseconds;
		timeToComplete = U.millisecondsToString(timeElapsed);
		milisecondsPerSpectrum = (double) timeElapsed / setSize;
		
		Peppy.Peppy.removeDuplicateMatches(positiveMatches);
		Peppy.Peppy.assignRankToMatches(positiveMatches);
		Match.setSortParameter(Match.SORT_BY_E_VALUE);
//		Match.setSortParameter(Match.SORT_BY_RANK_THEN_E_VALUE);
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
			if (match.getRank() == 1) {
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
		topReverseMatches = (new ScoringThreadServer(peptides, spectra, null)).getMatches();
		
		//Sort matches by e value	
		Match.setSortParameter(Match.SORT_BY_E_VALUE);
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
	
	/**
	 * This should be called only after our set of matches has been found
	 */
	private void calculateStastics() {
		//stats for tested matches
		boolean onePercenThresholdHasBeenReached = false;
		for (MatchContainer match: topForwardsTestedMatches) {
			if (match.isTrue()) {
				topRankTrueTally++;
			} else {
				topRankFalseTally++;
			}
			if (!onePercenThresholdHasBeenReached) {
				if ((double) topRankFalseTally / setSize >= 0.01) {
					onePercenThresholdHasBeenReached = true;
					trueTallyAtOnePercentError =  topRankTrueTally;
					percentAtOnePercentError = (double) trueTallyAtOnePercentError / setSize;
					eValueAtOnePercentError = match.getEValue();
				}
			}
		}
		//count total true
		for (MatchContainer match: testedMatches) {
			if (match.isTrue()) totalTrueTally++;
		}
		
		generatePrecisionRecallCurve();
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
//				recall = (double) trueCount / getSetSize();	
				recall = (double) trueCount / totalTrueTally;	
				
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
				File imageFile = new File(testDirectory, "/precision-recall.jpg");
				fileNameForPRCurve = testName + "/" + imageFile.getName();
				ImageIO.write(bdest,"JPG",imageFile);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

	private ArrayList<Peptide> loadCorrectPeptides() {
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
			correctAcidSequence = correctAcidSequence.trim();
			if (correctAcidSequence.equals("")) {
				validPeptideFile = false;
			}
			//adding to the array list
			if (validPeptideFile) {
				correctMatches.add(new Match(spectrum, new Peptide(correctAcidSequence), null));
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

	public int getTrueTallyAtOnePercentError() {
		return trueTallyAtOnePercentError;
	}

	public double getEValueAtOnePercentError() {
		return eValueAtOnePercentError;
	}

	public double getMilisecondsPerSpectrum() {
		return milisecondsPerSpectrum;
	}

	public double getAreaUnderPRCurve() {
		return areaUnderPRCurve;
	}

	public String getFileNameForPRCurve() {
		return fileNameForPRCurve;
	}

	public long getTimeElapsed() {
		return timeElapsed;
	}

	public String getTimeToComplete() {
		return timeToComplete;
	}

	public int getSetSize() {
		return setSize;
	}

	public double getPercentAtOnePercentError() {
		return percentAtOnePercentError;
	}

	public ArrayList<MatchContainer> getTestedMatches() {
		return testedMatches;
	}
	
	public double getEValueAtPercentForwards(double percent) {
		Match.setSortParameter(Match.SORT_BY_E_VALUE);
		Collections.sort(topForwardsMatches);
		int level = (int) (topForwardsMatches.size() * percent);
		return topForwardsMatches.get(level).getEValue();
	}
	
	public double getEValueAtPercentReverse(double percent) {
		Match.setSortParameter(Match.SORT_BY_E_VALUE);
		Collections.sort(topReverseMatches);
		int level = (int) (topReverseMatches.size() * percent);
		return topReverseMatches.get(level).getEValue();
	}


}
