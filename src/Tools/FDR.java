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
import Utilities.U;

public class FDR {
	
	public static void main(String args[]) {
		U.p("running FDR comparison report...");
		/* the tolerances we will be exploring */
		double [] precursorTolerances = {5, 10, 25};
		double [] fragmentTolerances = {10, 25, 50};
		
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
					area = findFDR(precursor, fragment, args);
					pw.println(precursor + "\t" + fragment + "\t" + area);
				}
			}
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		U.p("done");
	}
	
	/**
	 * 
	 * @param fragmentTolerance set to a negative number to ignore
	 * @param precursorTolerance set to a negative number to ignore
	 * @param args
	 */
	public static double findFDR(double precursorTolerance, double fragmentTolerance, String args[]) {
		
		U.startStopwatch();
		
		//grab our properties file, set up.
		Peppy.init(args);
		
		/* if we have specially defined tolerances, set those here */
		if (fragmentTolerance > 0) {
			Properties.fragmentTolerance = fragmentTolerance;
			U.p("fragmentTolerance: " + fragmentTolerance);
		}
		if (precursorTolerance > 0) {
			Properties.precursorTolerance = precursorTolerance;
			U.p("precursorTolerance: " + precursorTolerance);
		}
		NumberFormat nfPercent = NumberFormat.getPercentInstance();
		nfPercent.setMaximumFractionDigits(2);
		
		//Get references to our sequence files -- no nucleotide data is loaded at this point
		ArrayList<Sequence> sequences = Sequence_DNA.loadSequenceFiles(Properties.sequenceDirectoryOrFile);
		
		//This has to be 1 to properly calculate FDR
		Properties.maximumNumberOfMatchesForASpectrum = 1;
		
		//Loading a subset of our spectra
		U.p("loading spectral files...");
		ArrayList<File> spectraFiles = new ArrayList<File>();
		Spectrum.loadSpectraFilesFromFolder(Properties.spectraDirectoryOrFile, spectraFiles);
		U.p("loaded " + spectraFiles.size() + " spectra files");
		int setSize = Properties.numberOfSpectraToUseForFDR;
		if (setSize > spectraFiles.size()) setSize = spectraFiles.size();
		U.p("loading subset of spectra from the files...");
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
		
			
		U.p("running report for " + Properties.spectraDirectoryOrFile.getName());
		
		//getting forwards matches
		U.p("finding forwards matches...");
		ArrayList<Match> forwardsMatches = Peppy.getMatches(sequences, spectra);
		forwardsMatches = Peppy.reduceMatchesToOnePerSpectrum(forwardsMatches);
		Match.setSortParameter(Match.SORT_BY_E_VALUE);
		Collections.sort(forwardsMatches);
		
		//need to initialize things now that we have found matches
		sequences = Sequence_DNA.loadSequenceFiles(Properties.sequenceDirectoryOrFile);
		for (Spectrum spectrum: spectra) {
			spectrum.clearEValues();
		}
		
		//getting reverse matches -- need to reload the sequences
		U.p("finding reverse matches...");
		ArrayList<Match> reverseMatches = Peppy.getReverseMatches(sequences, spectra);
		reverseMatches = Peppy.reduceMatchesToOnePerSpectrum(reverseMatches);
		Match.setSortParameter(Match.SORT_BY_E_VALUE);
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
		Match.setSortParameter(Match.SORT_BY_E_VALUE);
		Collections.sort(allMatches);
		
		/* construct a list of points on the PR curve */
		ArrayList<Point2D.Double> points = new ArrayList<Point2D.Double>();
		int targetTally = 0;
		Match match;
		double precision, recall;
		for (int i = 0; i < allMatches.size(); i++) {
			match = allMatches.get(i);
			if (!match.isDecoy()) {
				targetTally++;
			}
			precision = (double) targetTally / (i + 1);
			recall = (double) targetTally / setSize;
			points.add(new Point2D.Double(precision, recall));
		}
		/* adjust for non-monotonic decreasing PR curve */
		double previousPrecision = points.get(points.size() - 1).x;
		for (int i = points.size() - 2; i >= 0; i--) {
			Point2D.Double point = points.get(i);
			if (point.x < previousPrecision) {
				point.x = previousPrecision;
			} else {
				previousPrecision = point.x;
			}
		}
		
		/* set a folder to store our reports */
		String identifierString = Properties.precursorTolerance + "-" + Properties.fragmentTolerance;
		File reportDir = new File("FDR/" + identifierString);
		reportDir.mkdirs();
		
		/* save curve image */
		PRCurve prCurve = new PRCurve(points);
		prCurve.writeFile(new File(reportDir, "pr-" + identifierString + ".jpg"));
		
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
			pw.println("Precision\tPercent Found\tE value");
			
			/* print the rest of the lines */
			double percent = 1;
			double increment = 0.01;
			for (int i = 1; i < points.size(); i++) {
				if (points.get(i).x < percent && points.get(i - 1).x != points.get(i).x) {
					percent -= increment;
					pw.println(points.get(i - 1).x + "\t" + points.get(i - 1).y + "\t" + allMatches.get(i - 1).getEValue());
				}
			}
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		U.stopStopwatch();
		
		return areaUndrePRCurve;
		
	}
	
	
	
	

}
