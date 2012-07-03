package Navigator;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

import Peppy.Properties;
import Peppy.U;

/**
 * A report that answers questions like:
 * 
 * what is the distribution of AAs?
 * what is the distribution of N-terminal AAs?
 * how do those differ from random distributions?
 * 
 * Copyright 2012, Brian Risk
 * 
 * @author Brian Risk
 *
 */
public class AAFrequency extends Report{
	
	
	public static void main(String args[]) {
		File reportFile = new File("/Users/risk2/Documents/workspace/JavaGFS/reports/WHIM16-Ellis033/6 WHIM16 - varimod/report.txt");
		File reportDirectory = new File("reports");
		AAFrequency report = new AAFrequency(loadMatches(reportFile));
		try {
			report.createReport(reportDirectory);
		} catch (IOException e) {
			e.printStackTrace();
		}
		U.p("done");
	}
	
	public AAFrequency(ArrayList<Match> matches) {
		super(matches);
	}
	

	@Override
	public File createReport(File reportDirectory) throws IOException {
		File acidFrequencyDirectory = new File(reportDirectory, "residue frequencies");
		acidFrequencyDirectory.mkdir();
		
		
		int totalAminoAcids = 0;
		int totalNTerminalAminoAcids = 0;
		
		/* get amino acid counts and store them in a hashtable */
		Hashtable<Character, Integer> acidsCounts = new Hashtable<Character, Integer>();
		Hashtable<Character, Integer> nTerminalAcidCounts = new Hashtable<Character, Integer>();
		for (Match match: matches) {
			String peptideSequence = match.getString("peptideSequence");
			char [] peptideAcids = peptideSequence.toCharArray();
			
			/* add in all acids except the n-terminal*/
			char acid;
			for (int i = 1; i < peptideAcids.length; i++) {
				acid = peptideAcids[i];
				/* ignore stops */
				if (acid == '.') continue;
				
				Integer acidCount = acidsCounts.get(acid);
				if (acidCount == null) {
					acidsCounts.put(acid, 1);
				} else {
					acidsCounts.put(acid, acidCount + 1);
				}	
				totalAminoAcids++;
			}
			
			/* add in just n-terminal acids */
			acid = peptideAcids[0];
			Integer acidCount = nTerminalAcidCounts.get(acid);
			if (acidCount == null) {
				nTerminalAcidCounts.put(acid, 1);
			} else {
				nTerminalAcidCounts.put(acid, acidCount + 1);
			}	
			totalNTerminalAminoAcids++;
			
		}
		
		
		File residueFrequencyFile = new File(acidFrequencyDirectory, "residue frequency.txt");
		PrintWriter residueFrequencyWriter = new PrintWriter(new FileWriter(residueFrequencyFile));
		ArrayList<Character> acidList = new ArrayList<Character>(acidsCounts.keySet());
		Collections.sort(acidList);
		for (Character acid: acidList) {
			int count = acidsCounts.get(acid);
			double percent = (double) count / totalAminoAcids;
			residueFrequencyWriter.println(acid + "\t" + count + "\t" + Properties.percentFormat.format(percent));
		}
		residueFrequencyWriter.flush();
		residueFrequencyWriter.close();
		
		
		File nTerminalFrequencyFile = new File(acidFrequencyDirectory, "n-terminal residue frequency.txt");
		PrintWriter nTerminalFrequencyWriter = new PrintWriter(new FileWriter(nTerminalFrequencyFile));
		acidList = new ArrayList<Character>(nTerminalAcidCounts.keySet());
		Collections.sort(acidList);
		for (Character acid: acidList) {
			int count = nTerminalAcidCounts.get(acid);
			double percent = (double) count / totalNTerminalAminoAcids;
			nTerminalFrequencyWriter.println(acid + "\t" + count + "\t" + Properties.percentFormat.format(percent) + "\t" + percent);
		}
		nTerminalFrequencyWriter.flush();
		nTerminalFrequencyWriter.close();
		
		
		File nTerminalFluctuationFile = new File(acidFrequencyDirectory, "n-terminal residue fluctuation.txt");
		PrintWriter nTerminalFluctuationWriter = new PrintWriter(new FileWriter(nTerminalFluctuationFile));
		acidList = new ArrayList<Character>(nTerminalAcidCounts.keySet());
		Collections.sort(acidList);
		for (Character acid: acidList) {
			int nTerminalCount = nTerminalAcidCounts.get(acid);
			double nTerminalPercent = (double) nTerminalCount / totalNTerminalAminoAcids;
			
			int totalCount = acidsCounts.get(acid);
			double totalPercent = (double) totalCount / totalAminoAcids;
			
			double percentDifference = nTerminalPercent - totalPercent;
			percentDifference *= 100;
			
			nTerminalFluctuationWriter.println(acid + "\t" + Properties.percentFormat.format(percentDifference));
		}
		nTerminalFluctuationWriter.flush();
		nTerminalFluctuationWriter.close();
		
		return acidFrequencyDirectory;
	}
	
	
	
	
	
	

}
