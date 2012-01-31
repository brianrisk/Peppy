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
import Peppy.Peppy;
import Peppy.Properties;
import Peppy.Sequence;
import Peppy.Sequence_DNA;
import Peppy.Spectrum;
import Peppy.U;

public class FDR {
	
	public static void main(String args[]) {
		U.p("running FDR comparison report...");
		U.startStopwatch();
		
		findFDR(-1, -1, args, Match.SORT_BY_IMP_VALUE);
//		iterate(args);
		
		U.stopStopwatch();
		U.p("done");
	}
	
	public static void getParameters() {
		
		
	}
	
	public static void iterate(String args[]) {
		/* precursor tolerances we will be exploring */
		U.p("List the precursor tolerances, separated by commas:");
		double [] precursorTolerances = extractArrayFromString(U.in());
		
		/* the fragment tolerances */
		U.p("List the fragment tolerances, separated by commas:");
		double [] fragmentTolerances = extractArrayFromString(U.in());
		
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
					area = findFDR(precursor, fragment, args, Match.SORT_BY_IMP_VALUE);
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
	public static double findFDR(double precursorTolerance, double fragmentTolerance, String args[], int sortParameter) {
		
		
		//grab our properties file, set up.
		Peppy.init(args);
		
		/* if we have specially defined tolerances, set those here */
		if (fragmentTolerance > 0) {
			Properties.fragmentTolerance = fragmentTolerance;
		}
		if (precursorTolerance > 0) {
			Properties.precursorTolerance = precursorTolerance;
		}
		NumberFormat nfPercent = NumberFormat.getPercentInstance();
		nfPercent.setMaximumFractionDigits(2);
		
		//Get references to our sequence files -- no nucleotide data is loaded at this point
		ArrayList<Sequence> sequences = Sequence_DNA.loadSequenceFiles(Properties.sequenceDirectoryOrFile);
		
		//This has to be 1 to properly calculate FDR
		Properties.maximumNumberOfMatchesForASpectrum = 1;
		
		//Loading a subset of our spectra
		ArrayList<File> spectraFiles = new ArrayList<File>();
		Spectrum.loadSpectraFilesFromFolder(Properties.spectraDirectoryOrFile, spectraFiles);
		int setSize = Properties.numberOfSpectraToUseForFDR;
		if (setSize > spectraFiles.size()) setSize = spectraFiles.size();
		ArrayList<Spectrum> spectra = new ArrayList<Spectrum>();
		Random random = new Random();
		File spectrumFile;
		while (spectra.size() < setSize) {
			spectrumFile = spectraFiles.remove(random.nextInt(spectraFiles.size()));
			spectra.addAll(Spectrum.loadSpectra(spectrumFile));
		}
		for (int i = 0; i < spectra.size(); i++) {
			spectra.get(i).setId(i);
		}
				
		//getting forwards matches
		ArrayList<Match> forwardsMatches = Peppy.getMatches(sequences, spectra);
		forwardsMatches = Peppy.reduceMatchesToOnePerSpectrum(forwardsMatches);
		Match.setSortParameter(sortParameter);
		Collections.sort(forwardsMatches);
		
		//need to initialize things now that we have found matches
		sequences = Sequence_DNA.loadSequenceFiles(Properties.sequenceDirectoryOrFile);
		for (Spectrum spectrum: spectra) {
			spectrum.clearEValues();
		}
		
		//getting reverse matches -- need to reload the sequences
		ArrayList<Match> reverseMatches = Peppy.getReverseMatches(sequences, spectra);
		reverseMatches = Peppy.reduceMatchesToOnePerSpectrum(reverseMatches);
		Match.setSortParameter(sortParameter);
		Collections.sort(reverseMatches);
		
		/* label the reverse matches as being from a reverse database */
		for (Match match: reverseMatches) {
			match.setDecoy(true);
		}
		
		/* combine the matches */
		ArrayList<Match> allMatches = new ArrayList<Match>(forwardsMatches.size() + reverseMatches.size());
		allMatches.addAll(forwardsMatches);
		allMatches.addAll(reverseMatches);
		Peppy.reduceMatchesToOnePerSpectrum(allMatches);
		
		/* sort the matches */
		Match.setSortParameter(sortParameter);
		Collections.sort(allMatches);
		
		/* construct a list of points on the PR curve */
		ArrayList<Point2D.Double> points = new ArrayList<Point2D.Double>();
		int totalCorrect = 0;
		double precision, recall;
		for (int i = 0; i < allMatches.size(); i++) {
			Match match = allMatches.get(i);
			if (!match.isDecoy()) {
				totalCorrect++;
			}
			precision = (double) totalCorrect / (i + 1);
			recall = (double) totalCorrect / setSize;
			points.add(new Point2D.Double(precision, recall));
		}
		
		
		/* set a folder to store our reports */
		String identifierString = Properties.precursorTolerance + "-" + Properties.fragmentTolerance;
		File reportDir = new File("FDR/" + identifierString);
		reportDir.mkdirs();
		
		/* save curve image */
		PRCurve prCurve = new PRCurve(points);
		prCurve.writeFile(new File(reportDir, "pr-" + identifierString + ".jpg"));
		
		/* print some basic stats */
		int totalComparisons = 0;
		for (Match match: allMatches)  {
			totalComparisons += match.getSpectrum().getEValueCalculator().getNumberOfMatches();
		}
		double correctPerComparison = (double) totalCorrect / totalComparisons;
		U.p(totalComparisons + "\t" + totalCorrect + "\t" + correctPerComparison);
		
		/* save FDR report */
		double areaUndrePRCurve = -1;
		try {
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
				if (points.get(i).x < percent && points.get(i - 1).x != points.get(i).x) {
					percent -= increment;
					pw.println(points.get(i - 1).x + "\t" + points.get(i - 1).y + "\t" + allMatches.get(i - 1).getIMP());
				}
			}
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
//		U.stopStopwatch();
		
		return areaUndrePRCurve;
		
	}
	
	private static void printMatches(ArrayList<Match> matches, String fileName) {
		File reportFile = new File(fileName);
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(reportFile)));

			for (int i = 0; i < matches.size(); i++) {
				Match match = matches.get(i);
				pw.println(match.toString());
			}

			pw.flush();
			pw.close();

		} catch (IOException e) {
			U.p("ERROR: Could not create file writer for our report");
			e.printStackTrace();
		}
	}
	
	
	
	

}
