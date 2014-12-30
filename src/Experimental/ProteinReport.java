package Experimental;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;

import Navigator.HTML;
import Navigator.Match;
import Navigator.MatchTable;
import Navigator.MatchesToProtein;
import Navigator.ProteinComparison;
import Peppy.U;


/**
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class ProteinReport {
	
	private double proteinScoreCutoff = 10;
//	private double matchScoreCutoff = 18.169118251665612;
//	private double matchScoreCutoff = 16.224111990961575;
//	private double matchScoreCutoff = 14;
	private double matchScoreCutoff = 22.556025416890137;
	
	
	
	
	public static void main(String args[]) {
		new ProteinReport();
	}
	
	
	public ProteinReport() {
//		mayo();
		WashU();
		U.p("done");
	}
	
	
	public void WashU() {
		Hashtable<String, MatchesToProtein> whim2 = new Hashtable<String, MatchesToProtein>();
		Hashtable<String, MatchesToProtein> whim16 = new Hashtable<String, MatchesToProtein>();
		
		/* 
		 * HACK to get protein locations
		 */
		ArrayList<Match> genomeMatches = Match.loadMatches(new File("/Users/risk2/PeppyData/WashU/reports/43/WHIM16/WashU WHIM16 xeno/report.txt"));
		genomeMatches.addAll(Match.loadMatches(new File("/Users/risk2/PeppyData/WashU/reports/43/WHIM16/WashU WHIM16 xeno/report.txt")));
		MatchTable matchTable = new MatchTable(true);
		for (Match match: genomeMatches) matchTable.put(match.getString("peptideSequence"), match);
		Hashtable<String, ArrayList<Match>> genomeHash = matchTable.getHashtable();
		
		/* 33 */
		addToProteinHash(new File("//Users/risk2/PeppyData/WashU/reports/33/WHIM2/WashU WHIM2 human/report.txt"), whim2);
		addToProteinHash(new File("//Users/risk2/PeppyData/WashU/reports/33/WHIM16/WashU WHIM16 human/report.txt"), whim16);
		
		/* 41 */
		addToProteinHash(new File("//Users/risk2/PeppyData/WashU/reports/41/WHIM2/WashU WHIM2 human/report.txt"), whim2);
		addToProteinHash(new File("//Users/risk2/PeppyData/WashU/reports/41/WHIM16/WashU WHIM16 human/report.txt"), whim16);
		
		/* 43 */
		addToProteinHash(new File("//Users/risk2/PeppyData/WashU/reports/43/WHIM2/WashU WHIM2 human/report.txt"), whim2);
		addToProteinHash(new File("//Users/risk2/PeppyData/WashU/reports/43/WHIM16/WashU WHIM16 human/report.txt"), whim16);
		
		
		/* link proteins to peptides */
		Hashtable<String, ArrayList<MatchesToProtein>> proteinsForEveryPeptideWHIM2 = getProteinsAssociatedWithPeptides(whim2);
		Hashtable<String, ArrayList<MatchesToProtein>> proteinsForEveryPeptideWHIM16 = getProteinsAssociatedWithPeptides(whim16);
		
		/* combine results */
		Hashtable<String, ProteinComparison> whims = new Hashtable<String, ProteinComparison>();
		Enumeration<MatchesToProtein> e;
		
		/* Integrating WHIM2 */
		e = whim2.elements();
		while (e.hasMoreElements()) {
			MatchesToProtein matchesToProtein = e.nextElement();
			ProteinComparison comparison = whims.get(matchesToProtein.getName());
			if (comparison == null) {
				comparison = new ProteinComparison(matchesToProtein.getName());
				whims.put(matchesToProtein.getName(), comparison);
			}
			comparison.addSet("whim2", matchesToProtein);
			
			/* determine overlaps with other proteins */
			matchesToProtein.determineProteinOverlaps(proteinsForEveryPeptideWHIM2);
		}
		
		/* Integrating WHIM16 */
		e = whim16.elements();
		while (e.hasMoreElements()) {
			MatchesToProtein matchesToProtein = e.nextElement();
			ProteinComparison comparison = whims.get(matchesToProtein.getName());
			if (comparison == null) {
				comparison = new ProteinComparison(matchesToProtein.getName());
				whims.put(matchesToProtein.getName(), comparison);
			}
			comparison.addSet("whim16", matchesToProtein);
			
			/* determine overlaps with other proteins */
			matchesToProtein.determineProteinOverlaps(proteinsForEveryPeptideWHIM16);
		}
		
		
		/* create a list and sort the list */
		ArrayList<ProteinComparison> comparisons = new ArrayList<ProteinComparison>(whims.values());
		Collections.sort(comparisons);
		
		/* remove comparisons from this list that are fully overlapped by other comparisons */
		for (int comparisonIndex = 0; comparisonIndex < comparisons.size(); comparisonIndex++) {
			if (comparisons.get(comparisonIndex).isFullyOverlapped()) {
				comparisons.remove(comparisonIndex);
				comparisonIndex--;
			}
		}
		
		/* save our reports */
		File reportDir = new File("protein comparison full");
		reportDir.mkdir();
		
		/* create main index file */
		try {

			PrintWriter matchWriter = new PrintWriter(new BufferedWriter(new FileWriter(new File(reportDir, "index.html"))));
			matchWriter.println(HTML.sortableTableHeader);
			matchWriter.println(HTML.tableTop);
			matchWriter.println("<tr>");
			matchWriter.println("<th>Uniprot</th>");
			matchWriter.println("<th>score</th>");
			matchWriter.println("<th>dominant</th>");
			matchWriter.println("<th>WHIM2</th>");
			matchWriter.println("<th>WHIM16</th>");
			matchWriter.println("<th>notes</th>");
			matchWriter.println("</tr>");
			matchWriter.println("</thead>");
			matchWriter.println("<tbody>");
			for (ProteinComparison comparison: comparisons) {
				
				/* skip comparisons hwere there is no unique peptide */
				if (comparison.getMaxUniquePeptideCount() == 0) continue;
				
//				/* So fucking hackey.  This tries to tell if the protein is in chr8 */
//				boolean in8 = false;
//				Protein dominant = comparison.getDominantProtein();
//				ArrayList<Match> dominantMatches = dominant.getMatches();
//				for (Match match: dominantMatches) {
//					ArrayList<Match> genomeHits = genomeHash.get(match.getString("peptideSequence"));
//					if (genomeHits == null) continue;
//					for (Match hitMatch: genomeHits) {
//						if (hitMatch.getInt("RankCount") == 1 && hitMatch.getString("sequenceName").equals("chr8")) {
//							in8 = true;
//							continue;
//						}
//					}
//				}
//				if (!in8) continue;
				
				
//					if (comparison.getMinMatchCount() < 20 && !comparison.hasUnique()) continue;
//				if (comparison.getMinMatchCount() < 10) continue;
				if (comparison.getScore() < 1.3) continue;
				String uniprot = "<a href=\"http://www.uniprot.org/uniprot/" + comparison.getName() + "\">";
				matchWriter.println("<tr>");
				matchWriter.println("<td>" + uniprot + comparison.getName() + "</a></td>");
				matchWriter.println("<td>" + comparison.getScore() + "</td>");
				matchWriter.println("<td>" + comparison.getDominantSampleName() + "</td>");
				
				MatchesToProtein matchesToProtein =  whim2.get(comparison.getName());
				if (matchesToProtein != null) {
					String fileName = "whim2-" + matchesToProtein.getName() + ".html";
					File saveFile = new File(reportDir, fileName);
					saveIndividualProteinReport(matchesToProtein, saveFile);
					
					String unique = "";
//						if (protein.hasUniqueMatch()) unique = " (" + protein.getUniquePeptideCount() + ")";
					matchWriter.println("<td><a href=\"" + fileName + "\">" + matchesToProtein.getMatches().size() + unique + "</a></td>");
//						matchWriter.println("<td>" + protein.getMatches().size() + unique + "</td>");
				} else {
					matchWriter.println("<td>0</td>");
				}
				
				
				matchesToProtein =  whim16.get(comparison.getName());
				if (matchesToProtein != null) {
					String fileName = "whim16-" + matchesToProtein.getName() + ".html";
					File saveFile = new File(reportDir, fileName);
					saveIndividualProteinReport(matchesToProtein, saveFile);
					
					String unique = "";
//						if (protein.hasUniqueMatch()) unique = " (" + protein.getUniquePeptideCount() + ")";
					matchWriter.println("<td><a href=\"" + fileName + "\">" + matchesToProtein.getMatches().size() + unique + "</a></td>");
//						matchWriter.println("<td>" + protein.getMatches().size() + unique + "</td>");
				} else {
					matchWriter.println("<td>0</td>");
				}
				

				
				matchWriter.println("<td></td>");
				matchWriter.println("</tr>");
			}
			matchWriter.println(HTML.sortableTableFooter);
			matchWriter.flush();
				matchWriter.close();
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
		
	}
	
	
	
	public void mayo() {
		U.p("loading results");

		File resultsFile = new File("/Users/risk2/PeppyData/Mayo/reports/human-mods/report.txt");
		
		Hashtable<String, MatchesToProtein> matchesToProteins = getProteinHash(resultsFile);
		
		/* put our protein hash into a list */
		ArrayList<MatchesToProtein> proteinList = new ArrayList<MatchesToProtein>(matchesToProteins.values());
		Collections.sort(proteinList);

		
		
		
		/* save our reports */
		File reportDir = new File("proteins");
		reportDir.mkdir();
		
		/* create main index file */
		try {

				PrintWriter matchWriter = new PrintWriter(new BufferedWriter(new FileWriter(new File(reportDir, "index.html"))));
				matchWriter.println(HTML.sortableTableHeader);
				matchWriter.println(HTML.tableTop);
				matchWriter.println("<tr>");
				matchWriter.println("<th>protein</th>");
				matchWriter.println("<th>score</th>");
				matchWriter.println("<th>has unique</th>");
				matchWriter.println("<th>Uniprot</th>");
				matchWriter.println("</tr>");
				matchWriter.println("</thead>");
				matchWriter.println("<tbody>");
				for (MatchesToProtein matchesToProtein: proteinList) {
					if (matchesToProtein.getScore() >= proteinScoreCutoff) {
						String link = "<a href=\"" + matchesToProtein.getName() + ".html\">";
						String uniprot = "<a href=\"http://www.uniprot.org/uniprot/" + matchesToProtein.getName() + "\">";
						matchWriter.println("<tr>");
						matchWriter.println("<td>" + link + matchesToProtein.getName() + "</a></td>");
						matchWriter.println("<td>" + matchesToProtein.getScore() + "</td>");
						matchWriter.println("<td>" + matchesToProtein.hasUniqueMatch() + "</td>");
						matchWriter.println("<td>" + uniprot + matchesToProtein.getName() + "</a></td>");
						matchWriter.println("</tr>");
					}
				}
				matchWriter.println("</table></body></html>");
				matchWriter.flush();
				matchWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	
	
	public Hashtable<String, MatchesToProtein> getProteinHash(File resultsFile) {
		Hashtable<String, MatchesToProtein> matchesToProteins = new Hashtable<String, MatchesToProtein>();
		addToProteinHash(resultsFile, matchesToProteins);
		return matchesToProteins;
		
	}
	
	public void addToProteinHash(File resultsFile,  Hashtable<String, MatchesToProtein> matchesToProteins) {
		/* load in our matches */
		ArrayList<Match> matches = Match.loadMatches(resultsFile);
		for (Match match: matches) { 
			match.set("reportFile", resultsFile);
		}
		for (int i = 0; i < matches.size(); i++) {
			if (matches.get(i).getScore() < matchScoreCutoff) {
				matches.remove(i);
				i--;
			}
		}
		
		/* we want to know how many hits each spectrum had.  
		 * This will inform us as to how much to "dilute" the score
		 * it would contribute to any given protein
		 */
		MatchTable matchesBySpectrum = new MatchTable(false);
		for (Match match: matches) {
			matchesBySpectrum.put(match.getString("spectrumMD5"), match);
		}
		
		/* add matches to the proteins, contribute to score */
		for (Match match: matches) {
			String proteinName = match.getString("SequenceName");
			MatchesToProtein existingProtein = matchesToProteins.get(proteinName);
			
			/* if that protein isn't in the list, make it */
			if (existingProtein == null) {
				existingProtein = new MatchesToProtein(proteinName);
				matchesToProteins.put(proteinName, existingProtein);
			}
			
			/* add this match to the protein */
			existingProtein.addMatch(match);
			
			/* see how many other places this particular peptide occurs */
			int numberOfPeptideOccurences = matchesBySpectrum.get(match.getString("spectrumMD5")).size();
			match.set("numberOfPeptideOccurences", numberOfPeptideOccurences);
			
			/*
			 * adding to the score.  We are normalizing for these two things:
			 * 
			 * 1) the number of times that peptide was found.
			 * 
			 * 2) the total match size.  This will allow for comparison between
			 * protein samples.  If sample one has 1 million spectra and sample two
			 * has 2 million, we don't want it to look like sample two has twice the
			 * expression of sample one;
			 */
//			existingProtein.addToScore((double) match.getString("peptideSequence").length() / (numberOfPeptideOccurences * matches.size()));
			existingProtein.addToScore(1.0 / (numberOfPeptideOccurences * matches.size()));
//			existingProtein.addToScore(1.0 / (matches.size()));
//			existingProtein.addToScore(1.0 / (numberOfPeptideOccurences));
			
		}
		
	}
	
	
	/**
	 * This returns a hash where the key is a peptide.  The values are all
	 * proteins that contain this peptide.
	 * 
	 * @param matchesToProteins this is a hash of all proteins stored with their names as key
	 * @return
	 */
	private Hashtable<String, ArrayList<MatchesToProtein>> getProteinsAssociatedWithPeptides(Hashtable<String, MatchesToProtein> matchesToProteins) {
		Hashtable<String, ArrayList<MatchesToProtein>> out = new Hashtable<String, ArrayList<MatchesToProtein>>();
		Enumeration<MatchesToProtein> e = matchesToProteins.elements();
		while (e.hasMoreElements()) {
			MatchesToProtein matchesToProtein = e.nextElement();
			ArrayList<Match> matches = matchesToProtein.getMatches();
			for (Match match: matches) {
				String peptide = match.getString("peptideSequence");
				ArrayList<MatchesToProtein> proteinsforPeptide = out.get(peptide);
				if (proteinsforPeptide == null) {
					proteinsforPeptide = new ArrayList<MatchesToProtein>();
					out.put(peptide, proteinsforPeptide);
				}
				proteinsforPeptide.add(matchesToProtein);
				
			}
		}
		return out;
	}
	
	
	private void saveIndividualProteinReport(MatchesToProtein matchesToProtein, File file) {
		
		
		try {
			PrintWriter matchWriter = new PrintWriter(new BufferedWriter(new FileWriter(file)));
			matchWriter.println(HTML.sortableTableHeader);
			
			matchWriter.println("<h1>" + matchesToProtein.getName() + "</h1>");
			
			matchWriter.println("<p>Number of identified peptides unique in proteome: " + matchesToProtein.getUniquePeptideCount() + "</p>");
			
			/* print the related proteins overlap table */
			ArrayList<ProteinOverlap> proteinOverlaps = matchesToProtein.getProteinOverlaps();
			matchWriter.println("<table class=\"tablesorter\"><thead>");
			matchWriter.println("<tr>");
			matchWriter.println("<th>protein</th>");
			matchWriter.println("<th>only in me</th>");
			matchWriter.println("<th>in both</th>");
			matchWriter.println("<th>only in other</th>");
			matchWriter.println("</tr>");
			matchWriter.println("</thead>");
			matchWriter.println("<tbody>");
			for (ProteinOverlap proteinOverlap: proteinOverlaps) {
				String uniprot = "<a href=\"http://www.uniprot.org/uniprot/" + proteinOverlap.getProteinB().getName() + "\">";
				matchWriter.println("<tr>");
				matchWriter.println("<td>" + uniprot + proteinOverlap.getProteinB().getName() + "</a></td>");
				matchWriter.println("<td>" + proteinOverlap.onlyACount + "</td>");
				matchWriter.println("<td>" + proteinOverlap.overlapCount + "</td>");
				matchWriter.println("<td>" + proteinOverlap.onlyBCount + "</td>");
				matchWriter.println("</tr>");
			}
			matchWriter.println("</tbody>");
			matchWriter.println("</table>");
			
			/* print all matches */
			matchWriter.println(HTML.tableTop);
			matchWriter.println("<tr>");
			matchWriter.println("<th>peptide</th>");
			matchWriter.println("<th>score</th>");
			matchWriter.println("<th>mod index</th>");
			matchWriter.println("<th>isUnique</th>");
			matchWriter.println("<th>mod mass</th>");
			matchWriter.println("</tr>");
			matchWriter.println("</thead>");
			matchWriter.println("<tbody>");
			ArrayList<Match> matches = matchesToProtein.getMatches();
			for (Match match: matches) {
				
				if (match.getString("SequenceName").equals(matchesToProtein.getName())) {
					matchWriter.println("<tr>");
					File spectrumPage = new File(match.getFile("reportFile").getParent(), "spectra/" + match.getInt("spectrumID") + ".html");
					String peptideSequence = match.getString("peptideSequence");
					matchWriter.println("<td><a href=\"" + spectrumPage.getAbsolutePath() + "\">" + peptideSequence + "</a></td>");
					matchWriter.println("<td>" + Math.round(match.getScore()) + "</td>");
					
					if (match.getBoolean("isModified")) {
						int modIndex = match.getInt("modIndex");
						char modifiedAcid = peptideSequence.charAt(modIndex);
						matchWriter.println("<td>" + (match.getInt("start") + modIndex + 1) + " (" + modifiedAcid + ")</td>");
					} else {
						matchWriter.println("<td></td>");
					}
					
					
					if ( match.getInt("RankCount") == 1) {
						matchWriter.println("<td>yes</td>");
					} else {
						matchWriter.println("<td></td>");
					}
					
					matchWriter.println("<td>" + match.getDouble("modMass") + "</td>");
					matchWriter.println("</tr>");
				}
								
			}
			matchWriter.println(HTML.sortableTableFooter);
			matchWriter.flush();
			matchWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
