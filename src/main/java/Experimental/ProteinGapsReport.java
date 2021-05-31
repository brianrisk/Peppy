package Experimental;

import Navigator.Match;
import Peppy.Protein;
import Peppy.ProteinCoverage;
import Peppy.SequenceProtein;
import Peppy.U;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

/**
 * Finds proteins that have the following characteristics:
 * 
 * 1) are uniquely identifiable by a peptide
 * 2) have suspiciously large, contiguous gaps in coverage
 * 
 * @author Brian Risk
 *
 */
public class ProteinGapsReport {
	
	
	
	public static void main(String [] args) {		
		/* loading matches */
		ArrayList<Match> matches = new ArrayList<Match>();
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/Documents/MD Anderson/Visit results/group1/1 group1 - UniProt_Human_2012_03.fasta/report.txt")));
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/Documents/MD Anderson/Visit results/group2/1 group2 - UniProt_Human_2012_03.fasta/report.txt")));
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/Documents/MD Anderson/Visit results/group3/1 group3 - UniProt_Human_2012_03.fasta/report.txt")));
		
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/Documents/MD Anderson/Visit results/group4/1 group4 - UniProt_Human_2012_03.fasta/report.txt")));
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/Documents/MD Anderson/Visit results/group5/1 group5 - UniProt_Human_2012_03.fasta/report.txt")));
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/Documents/MD Anderson/Visit results/group6/1 group6 - UniProt_Human_2012_03.fasta/report.txt")));
		
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/PeppyData/ENCODE/GM12878/reports/GM maternal/1 spectra uncompressed - gencode11.fasta/report.txt")));
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/PeppyData/CPTAC/reports/CompRef_Proteome_BI/1 CompRef_Proteome_BI - gencode12.fasta/report.txt")));
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/PeppyData/CPTAC/reports/CompRef_Proteome_JHUC_P5AB/1 CompRef_Proteome_JHUC_P5AB - gencode12.fasta/report.txt")));
		
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/PeppyData/CPTAC/reports/CompRef_Proteome_PNNL/1 CompRef_Proteome_PNNL - gencode12.fasta/report.txt")));
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/PeppyData/CPTAC/reports/WashU_045_046_P5/1 WashU_P5 - gencode12.fasta/report.txt")));
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/PeppyData/CPTAC/reports/WashU_045_046_P6/1 WashU_P6 - gencode12.fasta/report.txt")));
		
		
		/* loading proteins */
		
//		String proteinFileName = "/Users/risk2/PeppyData/public/sequences/protein/gencode12.fasta";
		String proteinFileName = "/Users/risk2/PeppyData/public/sequences/protein/UniProt_Human_2012_03.fasta";
		matches.addAll(Match.loadMatches(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/wei-yang/1 spectra - UniProt_Human_2012_03.fasta/report.txt")));
		
		
		
		
		/* all the fraction names */
//		ArrayList<String> fractionLabels = new ArrayList<String>();
//		fractionLabels.add("01to21");
//		fractionLabels.add("22to24");
//		fractionLabels.add("25to26");
//		fractionLabels.add("27to28");
//		fractionLabels.add("29to30");
//		fractionLabels.add("31to32");
//		fractionLabels.add("33to34");
//		fractionLabels.add("35to36");
//		fractionLabels.add("37to38");
//		fractionLabels.add("39to41");
//		fractionLabels.add("42to43");
//		fractionLabels.add("44to45");
//		fractionLabels.add("46to47");
//		fractionLabels.add("48to49");
//		fractionLabels.add("50to51");
//		fractionLabels.add("52to53");
//		fractionLabels.add("54to55");
//		fractionLabels.add("56to57");
//		fractionLabels.add("58to59");
//		fractionLabels.add("60to62");
//		fractionLabels.add("63to65");
//		fractionLabels.add("66to68");
//		fractionLabels.add("69to77");
//		fractionLabels.add("78to84");
//		
//		for (String fraction: fractionLabels) {
//			U.p("processing for: " + fraction);
//			new ProteinGapsReport(matches, fraction);
//		}
		
		
		new ProteinGapsReport(matches, null,proteinFileName);
		
		U.p("done");
		
	}
	
	public ProteinGapsReport(ArrayList<Match> matches, String fraction, String proteinFileName) {
		
		/* load proteins */
		File proteinFile = new File(proteinFileName);
		Hashtable<String, ProteinCoverage> proteins = new Hashtable<String,ProteinCoverage>();
		SequenceProtein proteinDatabase = new SequenceProtein(proteinFile);
		ArrayList<Protein> proteinsLoaded = proteinDatabase.getProteinsFromDatabase(false, false);	
		for (Protein protein: proteinsLoaded) {	
			proteins.put(protein.getName(), new ProteinCoverage(protein));
		}
		
		/* see which peptides are found in more than one protein form */
		Hashtable<String, Boolean> uniquePeptides = new Hashtable<String, Boolean>();
		String spectrumMD5;
		for (Match match: matches) {
			spectrumMD5 = match.getString("spectrumMD5");
			Boolean found = uniquePeptides.get(spectrumMD5);
			if (found == null) {
				uniquePeptides.put(spectrumMD5, true);
			} else {
				uniquePeptides.put(spectrumMD5, false);
			}
		}
		
		/* for each match, get it's protein, then catalog coverage */
		int start, stop;
		String proteinName;
		for (Match match: matches) {
			if (fraction != null) {
				if (match.getFile("FilePath").getName().indexOf(fraction) == -1) continue;
			}
			start = match.getInt("start");
			stop = match.getInt("stop");
			proteinName = match.getString("SequenceName");
			spectrumMD5 = match.getString("spectrumMD5");
			ProteinCoverage protein = proteins.get(proteinName);
			
			/* in some odd cases, the protein name is blank in the database */
			if (protein == null) continue;
			if (uniquePeptides.get(spectrumMD5)) {
				protein.setIdentifiable(true);
			}
			protein.addCoverage(match);
		}
		
		/* calculate scores */
		ArrayList<ProteinCoverage> proteinList = new ArrayList<ProteinCoverage>(proteins.values());
		for (ProteinCoverage protein: proteinList) {
			protein.calculateScore();
		}
		
		
		/* reduce to identifiable */
		ArrayList<ProteinCoverage> identifiableProteins = new ArrayList<ProteinCoverage>(proteins.size());
		for (ProteinCoverage protein: proteinList) {
			if (protein.isIdentifiable()) 
				identifiableProteins.add(protein);
		}
		
		/* quick report */
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("protein gaps report " + fraction + ".txt")));
		
			
			Collections.sort(identifiableProteins);
			
			int maxReport = 10;
			if (identifiableProteins.size() < maxReport) maxReport = identifiableProteins.size();
			for (int i = 0; i < maxReport; i++) {
				ProteinCoverage protein = identifiableProteins.get(i);
				pw.println(protein.getName());
				pw.println("hits: " + protein.getHits());
				pw.println("score: " + Math.round(protein.getScore()));
				boolean [] coverage = protein.getCoverage();
				StringBuffer print = new StringBuffer();
				for (int index = 0; index < coverage.length; index++) {
					if (index % 50 == 0 && index != 0) {
						pw.println(print);
						print = new StringBuffer();
					} else {
						if (index % 10 == 0 && index != 0) {
							print.append(" ");
						}
					}
					
					if (coverage[index]) {
						print.append(protein.getAcidString().charAt(index));
					} else {
						print.append("-");
					}
				}
				pw.println(print);
				pw.println();
			}
			
			pw.flush();
			pw.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		/* make report on one protein of interest */
//		try {
//			String geneOfInterest = "Q96AG4";
//			
//			ProteinCoverage protein = proteins.get(geneOfInterest);
//			if (protein != null) {
//				PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(geneOfInterest + "-" + fraction + ".txt")));
//				pw.println(protein.getName());
//				boolean [] coverage = protein.getCoverage();
//				StringBuffer print = new StringBuffer();
//				for (int index = 0; index < coverage.length; index++) {
//					if (index % 50 == 0 && index != 0) {
//						pw.println(print);
//						print = new StringBuffer();
//					} else {
//						if (index % 10 == 0 && index != 0) {
//							print.append(" ");
//						}
//					}
//					
//					if (coverage[index]) {
//						print.append(protein.getAcidString().charAt(index));
//					} else {
//						print.append("-");
//					}
//				}
//				pw.println(print);
//				pw.println();
//			
//				pw.flush();
//				pw.close();
//			}
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
			
		
	}
	
	public void generateReport(ArrayList<ProteinCoverage> identifiableProteins) {
		
	}
	
	
}
