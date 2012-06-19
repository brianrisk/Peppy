package Validate;

import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

import Graphs.PRCurve;
import Peppy.Match;
import Peppy.MatchesSpectrum;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.ScoringServer;
import Peppy.Spectrum;
import Peppy.SpectrumLoader;
import Peppy.U;

public class TestSet {
	
	private String testName;
	private ArrayList<Spectrum> spectra;
	private ArrayList<MatchesSpectrum> spectraMatches;
	private ArrayList<Match> positiveMatches = new ArrayList<Match>();
	private ArrayList<Match> topReverseMatches = null;
	private ArrayList<MatchContainer> testedMatches = null;
	private double averageNumberOfPeaksPerSpectrum;
	
	//statistics
	private int falseTally = 0;
	private int trueTally = 0;
	private int trueTallyAtFivePercentError = -1;
	private double percentAtFivePercentError = -1;
	
	//PR curve
	private double areaUnderPRCurve = 0;
	
	long timeElapsed = 0;
	double milisecondsPerSpectrum;
	
	public TestSet(String testDirectoryName, String testName, ArrayList<Match> matches, ArrayList<Spectrum> spectra) {
		this(testDirectoryName, testName);
		positiveMatches = matches;
		this.spectra = spectra;
		
		/* set up where we will hold all of the matches for our spectra */
		spectraMatches = new ArrayList<MatchesSpectrum>(spectra.size());
		for (Spectrum spectrum: spectra) {
			spectraMatches.add(new MatchesSpectrum(spectrum));
		}
	}
	

	public TestSet(String testDirectoryName, String testName) {
		this.testName = testName;
		
		//load spectra for this test
		spectra = SpectrumLoader.loadSpectraFromFolder(testDirectoryName + "/" + testName + "/spectra");
		
		/* set up where we will hold all of the matches for our spectra */
		spectraMatches = new ArrayList<MatchesSpectrum>(spectra.size());
		for (Spectrum spectrum: spectra) {
			spectraMatches.add(new MatchesSpectrum(spectrum));
		}
		
		/* spectrum statistics */
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
		ArrayList<Match> matches = Peppy.Peppy.getMatchesWithPeptides(peptides, spectraMatches);
		long stopTimeMilliseconds = System.currentTimeMillis();
		timeElapsed += stopTimeMilliseconds - startTimeMilliseconds;
		
		/* add the matches to the full list */
		if (matches != null) {
			positiveMatches.addAll(matches);
		}
		
		
	}
	

	public ArrayList<Match> getPositiveMatches() {
		return positiveMatches;
	}


	public void resetTest() {
		positiveMatches = new ArrayList<Match>();
		falseTally = 0;
		trueTally = 0;
		trueTallyAtFivePercentError = 0;
		percentAtFivePercentError = 0;
		areaUnderPRCurve = 0;
	}
	
	/**
	 * This should be called only after our set of matches has been found
	 */
	public void calculateStastics() {		
		/* finding how much time per spectrum */
		milisecondsPerSpectrum = (double) timeElapsed / spectra.size();


		Match.setSortParameter(Match.SORT_BY_SCORE);
		Collections.sort(positiveMatches);
		
		
		/* Identify if each match is true or false */
		testedMatches = new ArrayList<MatchContainer>(positiveMatches.size());
		for (Match match: positiveMatches) {
			MatchContainer testedMatch = new MatchContainer(match);
			testedMatches.add(testedMatch);
		}
		
		/* sort the tested matches */
		Collections.sort(testedMatches);
		
		
		/* reduce the tested matches to one per spectrum, keeping the one that is correct if there is a correct one 
		 * else keep the match with the best score */
		ArrayList<MatchContainer> reducedTestedMatches = new ArrayList<MatchContainer>(testedMatches.size());
		for (Spectrum spectrum: spectra) {
			MatchContainer toAdd = null;
			for (MatchContainer mc: testedMatches) {
				if (mc.getMatch().getSpectrum().getId() == spectrum.getId()) {
					if (toAdd == null) {
						toAdd = mc; 
					}
					if (toAdd.getMatch().getScore() < mc.getMatch().getScore()) {
						toAdd = mc;
					}
					if (toAdd.getMatch().getScore() == mc.getMatch().getScore() && mc.isTrue()) {
						toAdd = mc;
					}
				}
			}
			if (toAdd != null) {
				reducedTestedMatches.add(toAdd);
			}
		}
		testedMatches = reducedTestedMatches;
		
		/* ensure the testedMatches are considering the correct peptide */
		if (Properties.scoringMethodName.equals("Peppy.Match_IMP")){
			for (MatchContainer mc: testedMatches) {
				if (mc.getMatch().getScore() < mc.getTrueMatch().getScore()) {
					mc.getMatch().setPeptide(mc.getTrueMatch().getPeptide());
					mc.getMatch().recalculateIMP();
					mc.getMatch().calculateScore();
					mc.validate();
				}
			}
		}
		
		
		/* again sort the tested matches */
		Collections.sort(testedMatches);
	
		
		//count total true
		double highestIMP = 0;
		for (MatchContainer match: testedMatches) {
			if (match.isTrue()) {
				trueTally++;
				if (match.getMatch().getIMP() > highestIMP) highestIMP = match.getMatch().getIMP();
			} else {
				falseTally++;
			}
			if ((double) falseTally / (trueTally + falseTally) <= 0.05) {
				trueTallyAtFivePercentError =  trueTally;
				percentAtFivePercentError = (double) trueTallyAtFivePercentError / spectra.size();
			}
		}
		
		generatePrecisionRecallCurve();
	}

	
	public BufferedImage generatePrecisionRecallCurve() {
		
		int trueCount = 0;
		double precision = 0; 
		double recall = 0;

		ArrayList<Point2D.Double> points = new ArrayList<Point2D.Double> ();
		
		for(int i = 0; i < testedMatches.size(); i++) {
			MatchContainer match = testedMatches.get(i);
			if (match.isTrue()) {
				trueCount++;
			}
						
			precision = (double) trueCount / (i + 1);
			recall = (double) trueCount / getSetSize();				

			Point2D.Double point =  new Point2D.Double(recall, precision);
			points.add(point);
			
		}

		PRCurve prCurve = new PRCurve(points);
		areaUnderPRCurve = prCurve.calculateAreaUnderCurve();
		return prCurve.generatePrecisionRecallCurve();
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
		
		/* sorting the peptides is VERY important */
		Collections.sort(correctPeptides);
		return correctPeptides;
	}
	
	
	
	@SuppressWarnings("unused")
	private ArrayList<Match> loadCorrectMatches() {
		ArrayList<Match> correctMatches = new ArrayList<Match>();
		for(MatchesSpectrum matchesSpectrum: spectraMatches) {
			//find the file for the correct peptide
			File spectrumFile = matchesSpectrum.getSpectrum().getFile();
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
				correctMatches.add(Properties.matchConstructor.createMatch(matchesSpectrum, peptide));
			}
			
//			if (spectrum.getFile().getName().equals("T10707_Well_H13_1768.77_19185.mgf..pkl")) {
//				U.p ("Valid? " + validPeptideFile);
//			}
		}
		return correctMatches;
	}
	
	
	/**
	 * Here the peptide list is probably from a reverse database.
	 * We are wanting to see what matches we find to a database that does not
	 * have the correct answer.
	 * @param peptides
	 */
	public void findTrueMatchesInAFalseDatabase(ArrayList<Peptide> peptides) {
		
		//get the matches
		ScoringServer scoringServer = new ScoringServer(peptides, spectraMatches);
		scoringServer.findMatches();
		topReverseMatches = Peppy.Peppy.getMatchesFromSpectraMatches(spectraMatches);
		
		//Sort matches by e value	
		Match.setSortParameter(Match.SORT_BY_SCORE);
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


	public String getName() {
		return testName;
	}

	public int getTrueTally() {
		return trueTally;
	}


	public int getTrueTallyAtFivePercentError() {
		return trueTallyAtFivePercentError;
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
		return spectra.size();
	}

	public ArrayList<MatchContainer> getTestedMatches() {
		return testedMatches;
	}
	




}
