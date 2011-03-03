package Validate;

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

import Peppy.Match;
import Peppy.Peppy;
import Peppy.Properties;
import Peppy.Sequence;
import Peppy.Spectrum;
import Utilities.U;

public class FPR {
	
	
	public static void main(String args[]) {
		U.p("running FPR report...");
		U.startStopwatch();
		
		//grab our properties file, set up.
		Peppy.init(args);
		NumberFormat nfPercent = NumberFormat.getPercentInstance();
		nfPercent.setMaximumFractionDigits(2);
		
		//What scoring mechanism?
		String scoreName = Properties.scoringMethodName;
		U.p("running report for " + scoreName);
		
		//Get references to our sequence files -- no nucleotide data is loaded at this point
		ArrayList<Sequence> sequences = Sequence.loadSequences(Properties.sequenceDirectoryOrFile);
		
		//Loading a subset of our spectra
		U.p("loading spectral files...");
		ArrayList<File> spectraFiles = new ArrayList<File>();
		Spectrum.loadSpectraFilesFromFolder(Properties.spectraDirectoryOrFile, spectraFiles);
		U.p("loaded " + spectraFiles.size() + " spectra files");
		int setSize = 10000;
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
		Match.setSortParameter(Match.SORT_BY_E_VALUE);
		Collections.sort(forwardsMatches);
		
		//need to initialize things now that we have found matches
		sequences = Sequence.loadSequences(Properties.sequenceDirectoryOrFile);
		for (Spectrum spectrum: spectra) {
			spectrum.clearEValues();
		}
		
		//getting reverse matches -- need to reload the sequences
		U.p("finding reverse matches...");
		ArrayList<Match> reverseMatches = Peppy.getReverseMatches(sequences, spectra);
		Match.setSortParameter(Match.SORT_BY_E_VALUE);
		Collections.sort(reverseMatches);
		
		
		//printing the matches
		printMatches(forwardsMatches, "FPR-" + scoreName + "-forwards.txt");
		printMatches(reverseMatches, "FPR-" + scoreName + "-reverse.txt");
			
		
		//Save FPRs
		ArrayList<Point2D.Double> points = new ArrayList<Point2D.Double>();
		int forwardsIndex = 0;
		Point2D.Double point;
		int forwardsSize = forwardsMatches.size();
		for (int reverseIndex = 0; reverseIndex < reverseMatches.size(); reverseIndex++) {
			if (forwardsIndex == forwardsSize) break;
			if (forwardsIndex < 0) forwardsIndex = 0;
			while (forwardsMatches.get(forwardsIndex).getEValue() < reverseMatches.get(reverseIndex).getEValue()) {
				forwardsIndex++;
				if (forwardsIndex == forwardsSize) break;
			}
			if (forwardsIndex == forwardsSize) break;
			if (forwardsIndex < 0) break;
			forwardsIndex--;
			point = new Point2D.Double(((double) (reverseIndex + 1) / (forwardsIndex + 1)), reverseMatches.get(reverseIndex).getEValue());
			points.add(point);
		}
		
		double fpr01 = getFPR(points, 0.01);
		double fpr05 = getFPR(points, 0.05);
		
		//Find the percent of total spectra found at 1%
		int total01 = 0;
		double percent01 = 0;
		for (int i = 0; i < forwardsMatches.size(); i++) {
			if (forwardsMatches.get(i).getEValue() <= fpr01) {
				total01++;
			} else {
				break;
			}
		}
		percent01 = (double) total01 / setSize;
		
		//Find the percent of total spectra found at 5%
		int total05 = 0;
		double percent05 = 0;
		for (int i = 0; i < forwardsMatches.size(); i++) {
			if (forwardsMatches.get(i).getEValue() <= fpr05) {
				total05++;
			} else {
				break;
			}
		}
		percent05 = (double) total05/ setSize;
		
		//find peptide spectrum mass differences
		double averageMassDifference = 0;
		double averageAbsMassDifference = 0;
		double lowestMassDifference = Double.MAX_VALUE;
		double highestMassDifference = Double.MIN_VALUE;
		double massDifference;
		for (int i = 0; i < total01; i++) {
			massDifference = forwardsMatches.get(i).getSpectrum().getMass() - forwardsMatches.get(i).getPeptide().getMass();
			averageMassDifference += massDifference;
			averageAbsMassDifference += Math.abs(massDifference);
			if (lowestMassDifference > massDifference) lowestMassDifference = massDifference;
			if (highestMassDifference < massDifference) highestMassDifference = massDifference;
		}
		averageMassDifference /= total01;
		averageAbsMassDifference /= total01;
			
		File fprFile = new File("FPR-" + scoreName + ".txt");
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(fprFile)));
			
			//print header
			pw.println("scoring system used: " + scoreName);
			pw.println("number of spectra in our set: " + spectra.size());
			pw.println("precursor tolerance: " + Properties.spectrumToPeptideMassError);
			pw.println("peak tolerance: " + Properties.peakDifferenceThreshold);
			pw.println();
			
			//print results
			pw.println("database: " + Properties.sequenceDirectoryOrFile.getName());
			pw.println("1% FPR: " + fpr01);
			pw.println("percent found at 1% FPR: " + nfPercent.format(percent01));
			pw.println("5% FPR: " + fpr05);
			pw.println("percent found at 5% FPR: " + nfPercent.format(percent05));
			pw.println();
			
			//print peptide/spectrum mass differences
			pw.println("average mass difference: " + averageMassDifference);
			pw.println("average absolute value of mass difference: " + averageAbsMassDifference);
			pw.println("lowest mass difference: " + lowestMassDifference);
			pw.println("higest mass difference: " + highestMassDifference);
			pw.println();

			pw.flush();
			pw.close();
			
		} catch (IOException e) {
			U.p("ERROR: Could not create file writer for our report");
			e.printStackTrace();
		}
		U.stopStopwatch();
		U.p("done");
	}
	
	private static double getFPR(ArrayList<Point2D.Double> points, double percent) {
		double out = -1;
		for (Point2D.Double point: points) {
			if (point.getX() < percent) out = point.getY();
		}
		return out;
	}
	
	private static void printMatches(ArrayList<Match> matches, String fileName) {
		File reportFile = new File(fileName);
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(reportFile)));
			
			int upper = 1000;
			if (matches.size() < upper) upper = matches.size();
			for (int i = 0; i < upper; i++) {
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
