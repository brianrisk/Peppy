package Peppy;

import Graphs.PRCurve;

import java.awt.geom.Point2D;
import java.io.*;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;


/**
 * Before a full search is performed, we should know the score cutoff dictated by the FDR.
 * This class should be instantiated before full searches and used to calculate FDR thresholds.
 * The information this class collects can then be used for score thresholds.
 * 
 * It will also save a report containing a PR curve, the list of thresholds and the properties
 * used.
 * 
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class FDR {
	
	ArrayList<Point2D.Double> points = new ArrayList<Point2D.Double>();
	ArrayList<MatchesSpectrum> spectraMatches;
	boolean usedFullSetOfSpectra = false;
	
	/* set for a given FDR */
	double falseDiscoveryRate;
	double scoreThreshold;
	
	

	/**
	 * Constructor where we do not specify a set of peptides, therefore we will be
	 * using the set of sequences from the Properties
	 * @param fullSetOfSpectra
	 */
	public FDR(ArrayList<Spectrum> fullSetOfSpectra) {
		this(fullSetOfSpectra, null);
	}
	
	public FDR(ArrayList<Spectrum> fullSetOfSpectra, ArrayList<Peptide> peptides) {
		
		/* set the score cutoff to a reasonably low point.
		 * This line should at least be before we make our MatchesSpectrum objects */
		Properties.minimumScore = 7;
		
		//Loading a subset of our spectra
		int setSize = Properties.numberOfSpectraToUseForFDR;
		if (setSize < 100) setSize = fullSetOfSpectra.size();
		if (setSize > fullSetOfSpectra.size()) setSize = fullSetOfSpectra.size();
		if (setSize == fullSetOfSpectra.size()) usedFullSetOfSpectra = true;
		ArrayList<Spectrum> spectra = null;
		
		
		if (setSize == fullSetOfSpectra.size()) {
			spectra = fullSetOfSpectra;
		} else {
			Random random = new Random();
			spectra = new ArrayList<Spectrum>(fullSetOfSpectra);
			while (spectra.size() > setSize) {
				spectra.remove(random.nextInt(spectra.size()));
			}
		}
		
		/* set up where we will hold all of the matches for our spectra */
		spectraMatches = new ArrayList<MatchesSpectrum>(spectra.size());
		for (Spectrum spectrum: spectra) {
			MatchesSpectrum matchesSpectrum = new MatchesSpectrum(spectrum);
			
			/* keep only best matches */
			matchesSpectrum.setWhatToKeep(Matches.KEEP_ONLY_BEST_MATCHES);
			
			spectraMatches.add(matchesSpectrum);
		}
		
		/* if we do not have peptides, load in the sequences specified in properties */
		if (peptides == null) {
			
			/* Get references to our sequence files -- no nucleotide data is loaded at this point */
			ArrayList<Sequence> sequences = Sequence.loadSequenceFiles(Properties.sequenceDirectoryOrFile);
				
			/* getting forwards matches */
			Peppy.getMatches(sequences, spectraMatches);
	
			/* need to initialize sequences for second pass */
			sequences = Sequence.loadSequenceFiles(Properties.sequenceDirectoryOrFile);		
			
			/* getting reverse matches */
			Peppy.getDecoyMatches(sequences, spectraMatches);
			
		/* we have a group of peptides.  Make a reverse database from them */
		} else {
			/* getting target matches */
			Peppy.getMatchesWithPeptides(peptides, spectraMatches);
			
			/* create our decoy database */
			ArrayList<Peptide> decoyPeptides = getDecoyPeptidesFromPeptides(peptides);
			
			/* getting decoy matches */
			Peppy.getMatchesWithPeptides(decoyPeptides, spectraMatches);
			
		}
		
		/* sort the matches by score */
		Collections.sort(spectraMatches);
		
		/* construct a list of points on the PR curve */
		int truePositiveCount = 0;
		int falsePositiveCount = 0;
		int decoyAssignmentCount = 0;
		double precision, recall;
		for (int i = 0; i < spectraMatches.size(); i++) {
			
			/* if there are no matches for a spectrum, then we've reached the end of useful matches */
			if (spectraMatches.get(i).getMatches().size() == 0) break;
			
			/* get one match for a given set of spectrum matches */
			Match match = spectraMatches.get(i).getMatches().get(0);
			
			/* Tally target / decoy */
			if (match.getPeptide().isDecoy()) {
				decoyAssignmentCount++;
			}
			falsePositiveCount = 2 * decoyAssignmentCount;
//			falsePositiveCount = decoyAssignmentCount;
			truePositiveCount = (i + 1) - falsePositiveCount;
			
			precision = (double) truePositiveCount / (truePositiveCount + falsePositiveCount);
			recall = (double) truePositiveCount / setSize;
			Point2D.Double point = new Point2D.Double(recall, precision);
			points.add(point);
		}
		
		
	}
	
	
	public double getScoreThreshold(double falseDiscoveryRate) {
		int bestIndex = getCutoffIndex(falseDiscoveryRate);	
		
		if (bestIndex == -1) {
			return -1;
		} else {
			return spectraMatches.get(bestIndex).getScore();
		}
	}
	
	public int getCutoffIndex(double falseDiscoveryRate) {
		int bestIndex = -1;
		int truePositiveCount = 0;
		int falsePositiveCount = 0;
		int decoyAssignmentCount = 0;
		double precision;
		for (int matchIndex = 0; matchIndex < spectraMatches.size(); matchIndex++) {
			/* if there are no matches for a spectrum, then we've reached the end of useful matches */
			if (spectraMatches.get(matchIndex).getMatches().size() == 0) break;
			
			Match match = spectraMatches.get(matchIndex).getMatches().get(0);
			if (match.getPeptide().isDecoy()) {
				decoyAssignmentCount++;
			}
			falsePositiveCount = 2 * decoyAssignmentCount;
//			falsePositiveCount = decoyAssignmentCount;
			truePositiveCount = (matchIndex + 1) - falsePositiveCount;
			
			precision = (double) truePositiveCount / (truePositiveCount + falsePositiveCount);
			if ((1.0 - precision) <= falseDiscoveryRate) bestIndex = matchIndex;
		}	
		
		return bestIndex;
		
	}
	
	
	
	
	
	public void saveReport(File reportDirectory) {
		/* do nothing if no data */
		if (points.size() == 0) return;
		
		/* set a folder to store our reports */
		String identifierString = Properties.precursorTolerance + "-" + Properties.fragmentTolerance;
		File reportDir = new File(reportDirectory, "FDR/" + identifierString + " " + U.getDateString());
		reportDir.mkdirs();
		
		/* save our properties */
		Properties.generatePropertiesFile(reportDir);
		
		/* save curve image */
		PRCurve prCurve = new PRCurve(points);
		prCurve.writeFile(new File(reportDir, "pr-" + identifierString + ".jpg"));
		
		double areaUndrePRCurve = -1;
		try {
			
			/* REPORT FOR CALCULATED CUTOFF POINTS */
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(reportDir, "FDR-" + identifierString + ".txt"))));
			/* parameters */
			pw.println("sequenceDirectoryOrFile: " +Properties.sequenceDirectoryOrFile.getAbsolutePath());
			pw.println("spectraDirectoryOrFile: " +Properties.spectraDirectoryOrFile.getAbsolutePath());
			pw.println("precursorTolerance: " +Properties.precursorTolerance);
			pw.println("fragmentTolerance: " +Properties.fragmentTolerance);
			pw.println();
			
			/*print nice stats */
			areaUndrePRCurve = prCurve.calculateAreaUnderCurve();
			pw.println("Area under PR curve: " + areaUndrePRCurve);
			pw.println();
			
			/*print first line*/
			pw.println("Precision\tPercentFound\tNumberFound\tIMP value");
			
			/* print the rest of the lines */
//			double percent = 1;
//			double increment = 0.01;
//			for (int i = 1; i < points.size(); i++) {
//				if (points.get(i).y < percent && points.get(i - 1).x != points.get(i).x) {
//					percent -= increment;
//					String percentFoundString = Properties.percentFormat.format((points.get(i - 1).y)) + "%";
//					pw.println(percentFoundString + "\t" + points.get(i - 1).x + "\t" + spectraMatches.get(i - 1).getScore());
//				}
//			}
//			/* to print out the final one */
//			String percentFoundString = Properties.percentFormat.format((points.get(points.size()  - 1).y)) + "%";
//			pw.println(percentFoundString + "\t" + points.get(points.size() - 1).x + "\t" + spectraMatches.get(points.size() - 1).getScore());
			
			double fdrLevel = 0;
			double fdrIncrement = 0.01;
			for (int i = 0; i < 20; i++) {
				double scoreThreshold = getScoreThreshold(fdrLevel);
				double numberFound = 0;
				for (MatchesSpectrum spectrumMatches: spectraMatches) {
					if (spectrumMatches.getScore() >= scoreThreshold) {numberFound++;}
					else {break;}
				}
				double percentFound = numberFound / spectraMatches.size();
				NumberFormat round = NumberFormat.getInstance();
				round.setMaximumFractionDigits(2);
				String percentFoundString = Properties.percentFormat.format(percentFound);
				pw.println(round.format(fdrLevel) + "\t" + percentFoundString + "\t" + (int) numberFound + "\t" + scoreThreshold);
				fdrLevel += fdrIncrement;
			}
			
			pw.flush();
			pw.close();
			
			/* FULL LIST OF MATCHES WITH TARGET/DECOY ASSIGNMENTS */
			pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(reportDir, "target-decoy-" + identifierString + ".txt"))));
			for (MatchesSpectrum spectrumMatch: spectraMatches) {
				/* if there are no matches for a spectrum, then we've reached the end of useful matches */
				if (spectrumMatch.getMatches().size() == 0) break;
				
				/* print the match */
				Match match = spectrumMatch.getMatches().get(0);
				pw.println(!match.getPeptide().isDecoy() + "\t" + match.toString());
			}
			pw.flush();
			pw.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}


	public ArrayList<MatchesSpectrum> getSpectraMatches() {
		return spectraMatches;
	}


	public boolean usedFullSetOfSpectra() {
		return usedFullSetOfSpectra;
	}
	
	/**
	 * Takes a list of peptides and makes them decoys.  Decoys are formed thusly:
	 * 
	 * target:  ABCDEFGHIJK
	 * Decoy:   JIHGFEDCBAK
	 * 
	 * @param peptides
	 * @return
	 */
	public ArrayList<Peptide> getDecoyPeptidesFromPeptides(ArrayList<Peptide> peptides) {
		ArrayList<Peptide> decoyPeptides = new ArrayList<Peptide>(peptides.size());
		int lengthMinusOne;
		for (Peptide peptide: peptides) {
			String target = peptide.getAcidSequenceString();
			StringBuffer decoy = new StringBuffer();
			lengthMinusOne =  target.length() - 1;
			decoy.append(target.subSequence(0,lengthMinusOne));
			decoy.reverse();
			decoy.append(target.charAt(lengthMinusOne));
			Peptide decoyPeptide = new Peptide(decoy.toString());
			decoyPeptide.setDecoy(true);
			decoyPeptides.add(decoyPeptide);
		}
		
		return decoyPeptides;
	}


	
	
	
	

}
