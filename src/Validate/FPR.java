package Validate;

import java.awt.geom.Point2D;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import Peppy.Match;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.ProteinDigestion;
import Peppy.ScoringThreadServer;
import Peppy.Spectrum;
import Reports.HTMLReporter;
import Utilities.U;

public class FPR {
	
	static final int setSize = 10000;
	
	public static void main(String args[]) {
		U.p("running FPR report...");
		
		Properties.maximumNumberOfMatchesForASpectrum = 10;
		Properties.numberOfMissedCleavages = 2;
		
		//What scoring mechanism?
		boolean useTandemFit = true;
		String scoreName = "TandemFit";
		if (useTandemFit) {
			Properties.defaultScore = Properties.DEFAULT_SCORE_TANDEM_FIT;
			Properties.spectrumToPeptideMassError = 2.0;
			Properties.peakDifferenceThreshold = 0.3;
			Properties.peakDifferenceThreshold = 0.2;
		} else {
			scoreName = "HMM_Score";
			Properties.defaultScore = Properties.DEFAULT_SCORE_HMM;
			HMMScore.HMMClass.HmmSetUp();
			Properties.highIntensityCleaning = true;
			Properties.spectrumToPeptideMassError = 0.1;
			Properties.peakDifferenceThreshold = 0.5;
		}
		U.p("running report for " + scoreName);
		
		ArrayList<File> databaseFiles = new ArrayList<File>();
		databaseFiles.add(new File("/Users/risk2/Documents/sprot/encode-data/annotation_sets/uniprot_human_2010_08/uniprot_sprot_varsplic.fasta"));
		databaseFiles.add(new File("/Users/risk2/Documents/sprot/encode-data/annotation_sets/uniprot_human_2010_09/uniprot_sprot_human.fasta"));
//		databaseFiles.add(new File("uniprot_sprot_human_varsplic.fasta"));
//		databaseFiles.add(new File("uniprot_sprot_human.fasta"));
		
		U.p("loading spectral files...");
		ArrayList<File> spectraFiles = new ArrayList<File>();
//		Spectrum.loadSpectraFilesFromFolder(new File("spectra encode membrane/GO_mem_FASP_dta20100628"), spectraFiles);
		Spectrum.loadSpectraFilesFromFolder(new File("/Users/risk2/PeppyOverflow/spectra encode membrane/GO_mem_FASP_dta20100628"), spectraFiles);
		U.p("loaded " + spectraFiles.size() + " spectra files");
		
		U.p("loading subset of spectra...");
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
		
		
		File fprFile = new File("FPR-" + scoreName + ".txt");
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(fprFile)));
			
			//print header
			pw.println("scoring system used: " + scoreName);
			pw.println("number of spectra in our set: " + spectra.size());
			pw.println("precursor tolerance: " + Properties.spectrumToPeptideMassError);
			pw.println("peak tolerance: " + Properties.peakDifferenceThreshold);
			pw.println();
			
			//go through each database
			for (File databaseFile: databaseFiles) {
				//set up results file
				File forwardsFile = new File(U.getFileNameWithoutSuffix(databaseFile) + "-forwards.txt");
				PrintWriter forwardsResults = new PrintWriter(new BufferedWriter(new FileWriter(forwardsFile)));
				File reverseFile = new File(U.getFileNameWithoutSuffix(databaseFile) + "-reverse.txt");
				PrintWriter reverseResults = new PrintWriter(new BufferedWriter(new FileWriter(reverseFile)));
				
				U.p("running report for " + databaseFile.getName());
				//get peptides, matches, and process results
				ArrayList<Peptide> peptides = ProteinDigestion.getReversePeptidesFromFASTA(databaseFile);
				ArrayList<Match> reverseMatches = (new ScoringThreadServer(peptides, spectra, null)).getMatches();
				Peppy.Peppy.assignRankToMatches(reverseMatches);
				Peppy.Peppy.assignConfidenceValuesToMatches(reverseMatches);
				Match.setSortParameter(Match.SORT_BY_E_VALUE);
				Collections.sort(reverseMatches);
				
				U.p("now finding forwards matches...");
				peptides = ProteinDigestion.getPeptidesFromFASTA(databaseFile);
				ArrayList<Match> forwardsMatches = (new ScoringThreadServer(peptides, spectra, null)).getMatches();
				Peppy.Peppy.assignRankToMatches(forwardsMatches);
				Peppy.Peppy.assignConfidenceValuesToMatches(forwardsMatches);
				Match.setSortParameter(Match.SORT_BY_E_VALUE);
				Collections.sort(forwardsMatches);
				
//				//print results to file
//				String folderPrefix = scoreName + " " + U.getFileNameWithoutSuffix(databaseFile);
//				File reportDir = new File(Properties.reportDirectory, folderPrefix + "-forwards");
//				HTMLReporter forwardsReport = new HTMLReporter(forwardsMatches, spectra, null, reportDir);
//				forwardsReport.generateFullReport();
////				U.p("creating text reports");
////				TextReporter textReport = new TextReporter(forwardsMatches, spectra, null, reportDir);
////				textReport.generateFullReport();
//				
//				File reverseDir = new File(Properties.reportDirectory, folderPrefix + "-reverse");
//				HTMLReporter reverseReport = new HTMLReporter(reverseMatches, spectra, null, reverseDir);
//				reverseReport.generateFullReport();
////				U.p("creating text reports");
////				TextReporter reverseTextReport = new TextReporter(reverseMatches, spectra, null, reverseDir);
////				reverseTextReport.generateFullReport();
				
				
				for (Match match: forwardsMatches) {
					forwardsResults.println(match.getEValue() + "\t" + match.getPeptide().getAcidSequenceString());
				}
				forwardsResults.flush();
				forwardsResults.close();
				for (Match match: reverseMatches) {
					reverseResults.println(match.getEValue() + "\t" + match.getPeptide().getAcidSequenceString() + "\t" + match.getSpectrum().getFile().getName());
				}
				reverseResults.flush();
				reverseResults.close();
				
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
				
//				int limit = points.size();
////				if (limit > 100) limit = 100;
//				for (int i = 0; i < limit; i++) {
//					point = points.get(i);
//					pw.println(point.getX() + ", " + point.getY());
//				}
				
				double fpr01 = getFPR(points, 0.01);
				double fpr05 = getFPR(points, 0.05);
				
				double percent01 = 0;
				for (int i = 0; i < forwardsMatches.size(); i++) {
					if (forwardsMatches.get(i).getEValue() <= fpr01) {
						percent01 =  i;
					} else {
						break;
					}
				}
				percent01 /= setSize;
				
				double percent05 = 0;
				for (int i = 0; i < forwardsMatches.size(); i++) {
					if (forwardsMatches.get(i).getEValue() <= fpr05) {
						percent05 =  i;
					} else {
						break;
					}
				}
				percent05 /= setSize;
				
				//print results
				pw.println("database: " + databaseFile.getName());
				pw.println("reverse database size: " + peptides.size());
				pw.println("1% FPR: " + fpr01);
				pw.println("percent found at 1% FPR: " + percent01);
				pw.println("5% FPR: " + fpr05);
				pw.println("percent found at 5% FPR: " + percent05);
				pw.println();
			}
			
			pw.flush();
			pw.close();
			
		} catch (IOException e) {
			U.p("ERROR: Could not create file writer for our report");
			e.printStackTrace();
		}
		U.p("done");
	}
	
	private static double getFPR(ArrayList<Point2D.Double> points, double percent) {
		double out = -1;
		for (Point2D.Double point: points) {
			if (point.getX() < percent) out = point.getY();
		}
		return out;
	}
	
//	private static double getFPR(ArrayList<Match> matches, double percent) {
//		int level = (int) (matches.size() * percent);
//		return matches.get(level).getEValue();
//	}

}
