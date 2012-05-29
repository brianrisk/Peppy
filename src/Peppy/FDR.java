package Peppy;

import java.awt.geom.Point2D;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import Graphs.PRCurve;


/**
 * Before a full search is performed, we should know the score cutoff dictated by the FDR.
 * This class should be instantiated before full searches and used to calculate FDR thresholds.
 * The information this class collects can then be used for score thresholds.
 * 
 * It will also save a report containing a PR curve, the list of thresholds and the properties
 * used.
 * 
 * 
 * @author Brian Risk
 *
 */
public class FDR {
	
	ArrayList<Point2D.Double> points = new ArrayList<Point2D.Double>();
	ArrayList<Match> allMatches;
	ArrayList<MatchesSpectrum> spectraMatches;

	
	public FDR(ArrayList<Spectrum> fullSetOfSpectra) {
		
		/* Get references to our sequence files -- no nucleotide data is loaded at this point */
		ArrayList<Sequence> sequences = Sequence.loadSequenceFiles(Properties.sequenceDirectoryOrFile);
		
		
		//Loading a subset of our spectra
		int setSize = Properties.numberOfSpectraToUseForFDR;
		if (setSize < 100) setSize = fullSetOfSpectra.size();
		if (setSize > fullSetOfSpectra.size()) setSize = fullSetOfSpectra.size();
		ArrayList<Spectrum> spectra = null;
		
		/* set the score cutoff to a reasonably low point */
		Properties.minimumScore = 7;
		
		
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
			
			/* keep only one best Match.
			 * This saves memory and ensures the "competition" Elias and Gygi prescribed
			 */
			matchesSpectrum.setWhatToKeep(Matches.KEEP_ONLY_BEST_MATCHES);
			
			spectraMatches.add(matchesSpectrum);
		}
				
		/* getting forwards matches */
		Peppy.getMatches(sequences, spectraMatches);

		
		
		/* need to initialize sequences for second pass */
		sequences = Sequence.loadSequenceFiles(Properties.sequenceDirectoryOrFile);		
		
		//getting reverse matches -- need to reload the sequences
		Peppy.getDecoyMatches(sequences, spectraMatches);

		
		/* combine the matches, keeping only one top match per spectrum */
		allMatches = new ArrayList<Match>(spectraMatches.size());
		ArrayList<Match> matches;
		for (MatchesSpectrum matchesSpectrum: spectraMatches) {
			matches = matchesSpectrum.getMatches();
			if (matches.size() != 0) {
				allMatches.add(matches.get(0));
			}
		}
		
		/* sort the matches */
		Collections.sort(allMatches);
		
		/* construct a list of points on the PR curve */
		int truePositiveCount = 0;
		int falsePositiveCount = 0;
		int decoyAssignmentCount = 0;
		double precision, recall;
		for (int i = 0; i < allMatches.size(); i++) {
			Match match = allMatches.get(i);
			if (match.getPeptide().isDecoy()) {
				decoyAssignmentCount++;
			}
			falsePositiveCount = 2 * decoyAssignmentCount;
			truePositiveCount = (i + 1) - falsePositiveCount;
			
			precision = (double) truePositiveCount / (truePositiveCount + falsePositiveCount);
			recall = (double) truePositiveCount / setSize;
			Point2D.Double point = new Point2D.Double(recall, precision);
			points.add(point);
		}
		
		
	}
	
	
	public double getScoreThreshold(double falseDiscoveryRate) {
		int bestIndex = -1;
		int truePositiveCount = 0;
		int falsePositiveCount = 0;
		int decoyAssignmentCount = 0;
		double precision;
		for (int matchIndex = 0; matchIndex < allMatches.size(); matchIndex++) {
			Match match = allMatches.get(matchIndex);
			if (match.getPeptide().isDecoy()) {
				decoyAssignmentCount++;
			}
			falsePositiveCount = 2 * decoyAssignmentCount;
			truePositiveCount = (matchIndex + 1) - falsePositiveCount;
			
			precision = (double) truePositiveCount / (truePositiveCount + falsePositiveCount);
			if ((1.0 - precision) < falseDiscoveryRate) bestIndex = matchIndex;
		}	
		
		if (bestIndex == -1) {
			return -1;
		} else {
			return allMatches.get(bestIndex).getScore();
		}
	}
	
	
	
	public void saveReport(File reportDirectory) {
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
			pw.println("Precision\tPercent Found\tIMP value");
			
			/* print the rest of the lines */
			double percent = 1;
			double increment = 0.01;
			for (int i = 1; i < points.size(); i++) {
				if (points.get(i).y < percent && points.get(i - 1).x != points.get(i).x) {
					percent -= increment;
					String percentFoundString = Properties.nfPercent.format((points.get(i - 1).y)) + "%";
					pw.println(percentFoundString + "\t" + points.get(i - 1).x + "\t" + allMatches.get(i - 1).getScore());
				}
			}
			/* to print out the final one */
			String percentFoundString = Properties.nfPercent.format((points.get(points.size()  - 1).y)) + "%";
			pw.println(percentFoundString + "\t" + points.get(points.size() - 1).x + "\t" + allMatches.get(points.size() - 1).getScore());
			
			pw.flush();
			pw.close();
			
			/* FULL LIST OF MATCHES WITH TARGET/DECOY ASSIGNMENTS */
			pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(reportDir, "target-decoy-" + identifierString + ".txt"))));
			for (Match match: allMatches) {
				pw.println(!match.getPeptide().isDecoy() + "\t" + match.toString());
			}
			pw.flush();
			pw.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}
	

	
	
	
	

}
