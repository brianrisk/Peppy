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
import java.util.Hashtable;
import java.util.Random;

import Peppy.Match;
import Peppy.Peppy;
import Peppy.Properties;
import Peppy.Sequence;
import Peppy.Sequence_DNA;
import Peppy.Spectrum;
import Utilities.U;

public class FDR {
	
	public static final int Evalue = 1;
	public static final int IMP = 2;
	
	public static void main(String args[]) {
		U.p("running FDR report...");
		U.startStopwatch();
		
		//grab our properties file, set up.
		Peppy.init(args);
		NumberFormat nfPercent = NumberFormat.getPercentInstance();
		nfPercent.setMaximumFractionDigits(2);
		
		//What scoring mechanism?
		String scoreName = Properties.scoringMethodName;
		U.p("running report for " + scoreName);
		
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
		
		
		//printing the matches
		printMatches(forwardsMatches, "FDR-" + scoreName + "-forwards.txt");
		printMatches(reverseMatches, "FDR-" + scoreName + "-reverse.txt");
			
		
		//Save FPRs
		calculateFDR(IMP, forwardsMatches, reverseMatches, setSize);
		calculateFDR(Evalue, forwardsMatches, reverseMatches, setSize);
		
		
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
	
	private static void calculateFDR(int confidenceMethod, ArrayList<Match> forwardsMatches, ArrayList<Match> reverseMatches, int setSize) {
		//E value FPRs
		if (confidenceMethod == Evalue) Match.setSortParameter(Match.SORT_BY_E_VALUE);
		if (confidenceMethod == IMP) Match.setSortParameter(Match.SORT_BY_IMP_VALUE);
		Collections.sort(forwardsMatches);
		Collections.sort(reverseMatches);
		ArrayList<Point2D.Double> points = new ArrayList<Point2D.Double>();
		int forwardsIndex = 0;

		int forwardsSize = forwardsMatches.size();
		for (int reverseIndex = 0; reverseIndex < reverseMatches.size(); reverseIndex++) {
			if (forwardsIndex == forwardsSize) break;
			if (forwardsIndex < 0) forwardsIndex = 0;
			
			while (getConfidenceScore(forwardsMatches.get(forwardsIndex), confidenceMethod) < getConfidenceScore(reverseMatches.get(reverseIndex), confidenceMethod)) {
				forwardsIndex++;
				if (forwardsIndex == forwardsSize) break;
			}
			if (forwardsIndex == forwardsSize) break;
			if (forwardsIndex < 0) break;
			forwardsIndex--;
			Point2D.Double point = new Point2D.Double(((double) (reverseIndex + 1) / (forwardsIndex + 1)), getConfidenceScore(reverseMatches.get(reverseIndex), confidenceMethod));
			points.add(point);
		}
		
		double fpr01 = getFPR(points, 0.01);
		double fpr05 = getFPR(points, 0.05);
		
		//Find the percent of total spectra found at 1%
		int total01 = 0;
		double percent01 = 0;
		for (int i = 0; i < forwardsMatches.size(); i++) {
			if (getConfidenceScore(forwardsMatches.get(i), confidenceMethod) <= fpr01) {
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
			if (getConfidenceScore(forwardsMatches.get(i), confidenceMethod) <= fpr05) {
				total05++;
			} else {
				break;
			}
		}
		percent05 = (double) total05/ setSize;
		
		//find peptide spectrum mass differences
		double averageAbsMassDifference = 0;
		double lowestMassDifference = Double.MAX_VALUE;
		double highestMassDifference = Double.MIN_VALUE;
		double massDifference;
		for (int i = 0; i < total01; i++) {
			massDifference = forwardsMatches.get(i).getSpectrum().getMass() - forwardsMatches.get(i).getPeptide().getMass();
			averageAbsMassDifference += Math.abs(massDifference);
			if (lowestMassDifference > massDifference) lowestMassDifference = massDifference;
			if (highestMassDifference < massDifference) highestMassDifference = massDifference;
		}
		averageAbsMassDifference /= total01;
		
		//find number of spectra with matches in this range
		Hashtable<Integer, Integer> uniqueSpectrumIDsOnePercent = new Hashtable<Integer, Integer>();
		for (int i = 0; i < total01; i++) {
			Integer ID = new Integer(forwardsMatches.get(i).getSpectrum().getId());
			uniqueSpectrumIDsOnePercent.put(ID, ID);
		}
		Hashtable<Integer, Integer> uniqueSpectrumIDsFivePercent = new Hashtable<Integer, Integer>();
		for (int i = 0; i < total05; i++) {
			Integer ID = new Integer(forwardsMatches.get(i).getSpectrum().getId());
			uniqueSpectrumIDsFivePercent.put(ID, ID);
		}
			
		File fprFile = null;
		if (confidenceMethod == Evalue) fprFile = new File("FDR-E-" + Properties.scoringMethodName + ".txt");
		if (confidenceMethod == IMP) fprFile = new File("FDR-IMP-" + Properties.scoringMethodName + ".txt");
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(fprFile)));
			
			//set up percent formatting
			NumberFormat nfPercent = NumberFormat.getPercentInstance();
			nfPercent.setMaximumFractionDigits(2);
			
			//print header
			pw.println("scoring system used: " + Properties.scoringMethodName);
			pw.println("number of spectra in our set: " + setSize);
			pw.println("precursor tolerance: " + Properties.precursorTolerance);
			pw.println("peak tolerance: " + Properties.fragmentTolerance);
			pw.println();
			
			//print results
			pw.println("database: " + Properties.sequenceDirectoryOrFile.getName());
			pw.println("1% FDR: " + fpr01);
			pw.println("percent found at 1% FDR: " + nfPercent.format(percent01));
			pw.println("number of unique spectra at 1%: " + uniqueSpectrumIDsOnePercent.size());
			pw.println();
			pw.println("5% FDR: " + fpr05);
			pw.println("percent found at 5% FDR: " + nfPercent.format(percent05));
			pw.println("number of unique spectra at 5%: " + uniqueSpectrumIDsFivePercent.size());
			pw.println();
			
			//print peptide/spectrum mass differences
			pw.println("mass stats at 1%:");
			pw.println("average absolute value of mass difference: " + averageAbsMassDifference);
			pw.println("lowest mass difference: " + lowestMassDifference);
			pw.println("highest mass difference: " + highestMassDifference);
			pw.println();
			
			

			pw.flush();
			pw.close();
			
			/* saves a list of the masses along with the ion count */
//			if (confidenceMethod == Evalue) {
//				pw = new PrintWriter(new BufferedWriter(new FileWriter("ioncount.txt")));
//				for (Match match: forwardsMatches) {
//					if (match.getEValue() < fpr01) {
//						pw.println(match.getPeptide().getMass() + "," + match.getIonMatchTally());
//					}
//				}
//			}
//			pw.flush();
//			pw.close();
			
		} catch (IOException e) {
			U.p("ERROR: Could not create file writer for our report");
			e.printStackTrace();
		}
	}
	
	private static double getConfidenceScore(Match match, int confidenceMethod) {
		if (confidenceMethod == Evalue) return match.getEValue();
		if (confidenceMethod == IMP) return match.getIMP();
		return -1;
	}

}
