package Tools;

import java.awt.geom.Point2D;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import Graphs.PRCurve;
import Peppy.Match;
import Peppy.MatchesSpectrum;
import Peppy.Peppy;
import Peppy.Properties;
import Peppy.Sequence;
import Peppy.Spectrum;
import Peppy.SpectrumLoader;
import Peppy.U;

public class FDR {
	
	public static void main(String args[]) {
		
		U.startStopwatch();
		U.p("performing FDR analysis...");
		
//		if (args.length == 0) {
			findFDR(-1, -1, args);
//		} else {
//			iterate(args);
//		}
		
		U.stopStopwatch();
		U.p("done");
	}

	public FDR(ArrayList<Spectrum> spectra, ArrayList<Sequence> sequences) {
		
	}
	
	public static void iterate(String args[]) {
		/* precursor tolerances we will be exploring */
		U.p("List the precursor tolerances, separated by commas:");
		double [] precursorTolerances = extractArrayFromString(U.in());
		
		/* the fragment tolerances */
		U.p("List the fragment tolerances, separated by commas:");
		double [] fragmentTolerances = extractArrayFromString(U.in());
		
		U.p("processing...");
		/* the file where we will save the report summary */
		File reportDir = new File("FDR");
		reportDir.mkdir();
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(reportDir, "FDR comparison report.txt"))));
			pw.println("precursor\tfragment\tarea");
			double area, precursor, fragment;
			for (int precursorIndex = 0; precursorIndex < precursorTolerances.length; precursorIndex++) {
				for (int fragmentIndex = 0; fragmentIndex < fragmentTolerances.length; fragmentIndex++) {
					precursor = precursorTolerances[precursorIndex];
					fragment = fragmentTolerances[fragmentIndex];
					U.p("working on " + precursor + ", " + fragment);
					area = findFDR(precursor, fragment, args);
					pw.println(precursor + "\t" + fragment + "\t" + area);
				}
			}
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private static double [] extractArrayFromString (String string) {
		String [] split = string.split(",");
		double [] values = new double[split.length];
		for (int i = 0; i < split.length; i++) {
			values[i] = Double.parseDouble(split[i].trim());
		}
		return values;
	}
	
	/**
	 * 
	 * @param fragmentTolerance set to a negative number to ignore
	 * @param precursorTolerance set to a negative number to ignore
	 * @param args
	 */
	public static double findFDR(double precursorTolerance, double fragmentTolerance, String args[]) {
		
		
		//grab our properties file, set up.
		Peppy.init(args);
		
		/* if we have specially defined tolerances, set those here */
		if (fragmentTolerance > 0) {
			Properties.fragmentTolerance = fragmentTolerance;
		}
		if (precursorTolerance > 0) {
			Properties.precursorTolerance = precursorTolerance;
		}
		
		
		//Get references to our sequence files -- no nucleotide data is loaded at this point
		ArrayList<Sequence> sequences = Sequence.loadSequenceFiles(Properties.sequenceDirectoryOrFile);
		
		//Loading a subset of our spectra
		ArrayList<File> spectraFiles = new ArrayList<File>();
		SpectrumLoader.loadSpectraFilesFromFolder(Properties.spectraDirectoryOrFile, spectraFiles);
		int setSize = Properties.numberOfSpectraToUseForFDR;
		if (setSize > spectraFiles.size()) setSize = spectraFiles.size();
		ArrayList<Spectrum> spectra = new ArrayList<Spectrum>();
		Random random = new Random();
		File spectrumFile;
		while (spectra.size() < setSize) {
			spectrumFile = spectraFiles.remove(random.nextInt(spectraFiles.size()));
			spectra.addAll(SpectrumLoader.loadSpectra(spectrumFile));
		}
		
		/* set up where we will hold all of the matches for our spectra */
		ArrayList<MatchesSpectrum> spectraMatches = new ArrayList<MatchesSpectrum>(spectra.size());
		for (Spectrum spectrum: spectra) {
			spectraMatches.add(new MatchesSpectrum(spectrum));
		}
				
		//getting forwards matches
		U.p("getting forwards matches");
		ArrayList<Match> forwardsMatches = Peppy.getMatches(sequences, spectraMatches);
		forwardsMatches = Peppy.reduceMatchesToOnePerSpectrum(forwardsMatches);
		Collections.sort(forwardsMatches);
		
		//need to initialize things now that we have found matches
		sequences = Sequence.loadSequenceFiles(Properties.sequenceDirectoryOrFile);
		for (MatchesSpectrum matchesSpectrum: spectraMatches) {
			matchesSpectrum.clearMatches();
		}
		
		//getting reverse matches -- need to reload the sequences
		U.p("getting reverse matches");
		ArrayList<Match> reverseMatches = Peppy.getDecoyMatches(sequences, spectraMatches);
		reverseMatches = Peppy.reduceMatchesToOnePerSpectrum(reverseMatches);
		Collections.sort(reverseMatches);
		
		
		/* combine the matches */
		U.p("combining sets and analysing...");
		ArrayList<Match> allMatches = new ArrayList<Match>(forwardsMatches.size() + reverseMatches.size());
		allMatches.addAll(forwardsMatches);
		allMatches.addAll(reverseMatches);
		allMatches = Peppy.reduceMatchesToOnePerSpectrum(allMatches);
		
		/* sort the matches */
		Collections.sort(allMatches);
		
		/* construct a list of points on the PR curve */
		ArrayList<Point2D.Double> points = new ArrayList<Point2D.Double>();
		int truePositiveCount = 0;
		int falsePositiveCount = 0;
		int decoyAssignmentCount = 0;
		int targetAssignmentCount = 0;
		double precision, recall;
		for (int i = 0; i < allMatches.size(); i++) {
			Match match = allMatches.get(i);
			if (match.getPeptide().isDecoy()) {
				decoyAssignmentCount++;
			} else {
				targetAssignmentCount++;
			}
			falsePositiveCount = 2 * decoyAssignmentCount;
			truePositiveCount = (i + 1) - falsePositiveCount;
			
			precision = (double) truePositiveCount / (truePositiveCount + falsePositiveCount);
			recall = (double) truePositiveCount / setSize;
			Point2D.Double point = new Point2D.Double(recall, precision);
			points.add(point);
		}
		
		U.p();
		U.p("falsePositiveCount: " + falsePositiveCount);
		U.p("truePositiveCount: " + truePositiveCount);
		U.p("decoyAssignmentCount: " + decoyAssignmentCount);
		U.p("targetAssignmentCount: " + targetAssignmentCount);
		U.p("allMatches.size(): " + allMatches.size());
		U.p("setSize: " + setSize);
		
//		int totalCorrect = 0;
//		double precision, recall;
//		for (int i = 0; i < allMatches.size(); i++) {
//			Match match = allMatches.get(i);
//			if (!match.isDecoy()) totalCorrect++;
//			precision = (double) totalCorrect / (i + 1);
//			recall = (double) totalCorrect / setSize;
//			points.add(new Point2D.Double(precision, recall));
//		}
		
		
		/* set a folder to store our reports */
		String identifierString = Properties.precursorTolerance + "-" + Properties.fragmentTolerance;
		File reportDir = new File("FDR/" + identifierString + " " + U.getDateString());
		reportDir.mkdirs();
		
		/* save our properties */
		Properties.generatePropertiesFile(reportDir);
		
		/* save curve image */
		PRCurve prCurve = new PRCurve(points);
		prCurve.writeFile(new File(reportDir, "pr-" + identifierString + ".jpg"));
		
		/* save FDR report */
		double areaUndrePRCurve = -1;
		try {
			
			/* how we format our percents */
			NumberFormat nfPercent = NumberFormat.getPercentInstance();
			nfPercent.setMaximumFractionDigits(2);
			
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
					String percentFoundString = nfPercent.format((points.get(i - 1).y * 100)) + "%";
					pw.println(percentFoundString + "\t" + points.get(i - 1).x + "\t" + allMatches.get(i - 1).getScore());
				}
			}
			/* to print out the final one */
			String percentFoundString = nfPercent.format((points.get(points.size()  - 1).y * 100)) + "%";
			pw.println(percentFoundString + "\t" + points.get(points.size() - 1).x + "\t" + allMatches.get(points.size() - 1).getScore());
			
			pw.flush();
			pw.close();
			
			/* FULL LIST OF MATCHES WITH TARGET/DECOY ASSIGNMENTS */
			pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(reportDir, "target-decoy-" + identifierString + ".txt"))));
			for (Match match: allMatches) {
				pw.println(match.getPeptide().isDecoy() + "\t" + match.toString());
			}
			pw.flush();
			pw.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
//		U.stopStopwatch();
		
		return areaUndrePRCurve;
		
	}
	

	
	
	
	

}
