package Navigator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;

import Peppy.U;

/**
 * Copyright 2012, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class PeptideReport {
	
	private double matchScoreCutoff = 22.556025416890137;
	
	
	public static void main(String args[]) {
		new PeptideReport();
	}
	
	
	public PeptideReport() {
		U.p("creating report");
		createReport();
		U.p("done");
	}
	
	
	public void createReport() {
		Hashtable<String, MatchesToPeptide> whim2 = new Hashtable<String, MatchesToPeptide>();
		Hashtable<String, MatchesToPeptide> whim16 = new Hashtable<String, MatchesToPeptide>();
		
		/* 
		 * HACK to get peptide locations
		 */
//		ArrayList<Match> genomeMatches = Match.loadMatches(new File("/Users/risk2/PeppyData/WashU/reports/43/WHIM16/WashU WHIM16 xeno/report.txt"));
//		genomeMatches.addAll(Match.loadMatches(new File("/Users/risk2/PeppyData/WashU/reports/43/WHIM16/WashU WHIM16 xeno/report.txt")));
//		MatchTable matchTable = new MatchTable(true);
//		for (Match match: genomeMatches) matchTable.put(match.getString("peptideSequence"), match);
//		Hashtable<String, ArrayList<Match>> genomeHash = matchTable.getHashtable();
		
		/* 33 */
		addToHash(new File("//Users/risk2/PeppyData/WashU/reports/33/WHIM2/WashU WHIM2 hg19/report.txt"), whim2);
		addToHash(new File("//Users/risk2/PeppyData/WashU/reports/33/WHIM16/WashU WHIM16 hg19/report.txt"), whim16);
		
		/* 41 */
		addToHash(new File("//Users/risk2/PeppyData/WashU/reports/41/WHIM2/WashU WHIM2 hg19/report.txt"), whim2);
		addToHash(new File("//Users/risk2/PeppyData/WashU/reports/41/WHIM16/WashU WHIM16 hg19/report.txt"), whim16);
		
		/* 43 */
		addToHash(new File("//Users/risk2/PeppyData/WashU/reports/43/WHIM2/WashU WHIM2 hg19/report.txt"), whim2);
		addToHash(new File("//Users/risk2/PeppyData/WashU/reports/43/WHIM16/WashU WHIM16 hg19/report.txt"), whim16);
		
		
		
		/* combine results */
		Hashtable<String, PeptideComparison> whims = new Hashtable<String, PeptideComparison>();
		Enumeration<MatchesToPeptide> e;
		
		/* Integrating WHIM2 */
		e = whim2.elements();
		while (e.hasMoreElements()) {
			MatchesToPeptide matchesToPeptide = e.nextElement();
			PeptideComparison comparison = whims.get(matchesToPeptide.getName());
			if (comparison == null) {
				comparison = new PeptideComparison(matchesToPeptide.getName());
				whims.put(matchesToPeptide.getName(), comparison);
			}
			comparison.addSet("whim2", matchesToPeptide);
			
		}
		
		/* Integrating WHIM16 */
		e = whim16.elements();
		while (e.hasMoreElements()) {
			MatchesToPeptide matchesToPeptide = e.nextElement();
			PeptideComparison comparison = whims.get(matchesToPeptide.getName());
			if (comparison == null) {
				comparison = new PeptideComparison(matchesToPeptide.getName());
				whims.put(matchesToPeptide.getName(), comparison);
			}
			comparison.addSet("whim16", matchesToPeptide);
			
		}
		
		
		/* create a list and sort the list */
		ArrayList<SetComparison> comparisons = new ArrayList<SetComparison>(whims.values());
		Collections.sort(comparisons);
	
		
		/* save our reports */
		File reportDir = new File("peptide comparison");
		reportDir.mkdir();
		
		/* create main index file */
		try {

			PrintWriter matchWriter = new PrintWriter(new BufferedWriter(new FileWriter(new File(reportDir, "index.html"))));
			matchWriter.println(HTML.sortableTableHeader);
			matchWriter.println(HTML.tableTop);
			matchWriter.println("<tr>");
			matchWriter.println("<th>peptide</th>");
			matchWriter.println("<th>score</th>");
			matchWriter.println("<th>dominant</th>");
			matchWriter.println("<th>WHIM2</th>");
			matchWriter.println("<th>WHIM16</th>");
			matchWriter.println("<th>Uniprot</th>");
			matchWriter.println("<th>notes</th>");
			matchWriter.println("</tr>");
			matchWriter.println("</thead>");
			matchWriter.println("<tbody>");
			for (SetComparison comparison: comparisons) {
				
				if (comparison.getMatchTotal() < 200) continue;
				
//				if (comparison.getMinMatchCount() < 10) continue;
				if (comparison.getScore() < 1.3) continue;
				matchWriter.println("<tr>");
				matchWriter.println("<td>"  + comparison.getName() + "</td>");
				matchWriter.println("<td>" + comparison.getScore() + "</td>");
				matchWriter.println("<td>" + comparison.getDominantSampleName() + "</td>");
				
				MatchesToPeptide matchesToPeptide =  whim2.get(comparison.getName());
				if (matchesToPeptide != null) {
//						if (peptide.hasUniqueMatch()) unique = " (" + peptide.getUniquePeptideCount() + ")";
					matchWriter.println("<td>" + matchesToPeptide.getMatches().size() + "</td>");
//						matchWriter.println("<td>" + peptide.getMatches().size()  + "</td>");
				} else {
					matchWriter.println("<td>0</td>");
				}
				
				
				matchesToPeptide =  whim16.get(comparison.getName());
				if (matchesToPeptide != null) {
//						if (peptide.hasUniqueMatch()) unique = " (" + peptide.getUniquePeptideCount() + ")";
					matchWriter.println("<td>" + matchesToPeptide.getMatches().size() + "</td>");
//						matchWriter.println("<td>" + peptide.getMatches().size() + unique + "</td>");
				} else {
					matchWriter.println("<td>0</td>");
				}
				

				matchWriter.println("<td>" + comparison.getName() + "</td>");
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
		
	
	
	public Hashtable<String, MatchesToPeptide> getPeptideHash(File resultsFile) {
		Hashtable<String, MatchesToPeptide> matchesToPeptides = new Hashtable<String, MatchesToPeptide>();
		addToHash(resultsFile, matchesToPeptides);
		return matchesToPeptides;
	}
	
	
	public void addToHash(File resultsFile,  Hashtable<String, MatchesToPeptide> matchesToPeptides) {
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
		 * it would contribute to any given peptide
		 */
		MatchTable matchesBySpectrum = new MatchTable(false);
		for (Match match: matches) {
			matchesBySpectrum.put(match.getString("spectrumMD5"), match);
		}
		
		/* add matches to the peptides, contribute to score */
		for (Match match: matches) {
			String peptideName = match.getString("peptideSequence");
			MatchesToPeptide existingPeptide = matchesToPeptides.get(peptideName);
			
			/* if that peptide isn't in the list, make it */
			if (existingPeptide == null) {
				existingPeptide = new MatchesToPeptide(peptideName);
				matchesToPeptides.put(peptideName, existingPeptide);
			}
			
			/* add this match to the peptide */
			existingPeptide.addMatch(match);
			
			/* see how many other places this particular peptide occurs */
			int numberOfPeptideOccurences = matchesBySpectrum.get(match.getString("spectrumMD5")).size();
			match.set("numberOfPeptideOccurences", numberOfPeptideOccurences);
			
			/*
			 * adding to the score.  We are normalizing for these two things:
			 * 
			 * 1) the number of times that peptide was found.
			 * 
			 * 2) the total match size.  This will allow for comparison between
			 * peptide samples.  If sample one has 1 million spectra and sample two
			 * has 2 million, we don't want it to look like sample two has twice the
			 * expression of sample one;
			 */
			existingPeptide.addToScore(1.0 / (numberOfPeptideOccurences * matches.size()));

			
		}
		
	}
	

	
	
	

}
