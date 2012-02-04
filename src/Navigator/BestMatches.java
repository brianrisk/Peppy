package Navigator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;

import Peppy.U;

public class BestMatches {
	
	/* How results files are laid out
	 * 
	 *  0 ID
		1 SpectrumID
		2 MD5
		3 FileName
		4 Score
		5 PrecursorM/Z
		6 PrecursorNeutralMass
		7 E Value
		8 PeptideSequence
	 */
	
	/* our results files */
	String nucleotideFileName;
	String proteinFileName;
	
	/* where we keep the best results */
	Hashtable<String, Match_fromResults> bestResults = new Hashtable<String, Match_fromResults>();
	
	public static void main(String args[]) {
		BestMatches bestMatches = new BestMatches(
				"/Users/risk2/PeppyOverflow/WashU/reports/WHIM2 Xeno/spectra_1324954973534_report.txt",
				"/Users/risk2/PeppyOverflow/WashU/reports/WHIM2 Uniprot/WHIM2 unfrac and frac tremble/spectra_1324090184795_report.txt");
	}
	
	
	public BestMatches(String nucleotideFileName, String proteinFileName) {
		super();
		this.nucleotideFileName = nucleotideFileName;
		this.proteinFileName = proteinFileName;
		
		U.p("loading protein results");
		loadResults(proteinFileName, Match_fromResults.TYPE_PROTEIN);
		U.p("loading gnome results");
		loadResults(nucleotideFileName, Match_fromResults.TYPE_GENOME);
		
		/* reduce the best results down to the best one match for any given peptide */
		U.p("reducing and saving");
		Hashtable<String, Match_fromResults> bestPeptides = new Hashtable<String, Match_fromResults>();
		Enumeration<Match_fromResults> values = bestResults.elements();
		while (values.hasMoreElements()) {
			Match_fromResults match = values.nextElement();
			Match_fromResults otherMatch = bestPeptides.get(match.getPeptide());
			if (otherMatch == null) {
				bestPeptides.put(match.getPeptide(), match);
			} else {
				if (match.getScore() > otherMatch.getScore()) {
					bestPeptides.put(match.getPeptide(), match);
				}
			}
		}
		
		ArrayList<Match_fromResults> bestArray = new ArrayList<Match_fromResults>(bestPeptides.values());
		Collections.sort(bestArray);
		
		/* save all matches */
		try {
			PrintWriter allMatches = new PrintWriter(new BufferedWriter(new FileWriter("best matches.txt")));
			for (Match_fromResults match: bestArray) {
				allMatches.println(match.toString());
			}
			allMatches.flush();
			allMatches.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		/* save the matches unique to protein */
		try {
			PrintWriter bestGenome = new PrintWriter(new BufferedWriter(new FileWriter("genome unique matches.txt")));
			for (Match_fromResults match: bestArray) {
				if (match.getMatchType() == Match_fromResults.TYPE_GENOME) {
					bestGenome.println(match.toString());
				}
			}
			bestGenome.flush();
			bestGenome.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		U.p("done");
	}



	private void loadResults(String resultsFilePath, int resultsType) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(resultsFilePath)));
			
			/* read the header lines */
			br.readLine();
			br.readLine();
			br.readLine();
			br.readLine();
			
			/* read in the first line */
			String line = br.readLine();
			
			while (line != null) {
				String [] chunks = line.split("\t");
				
				/*
				 * (String peptide, String spectrumMD5, String spectrumFile,
			double score, int matchType) {
				 */
				Match_fromResults match = new Match_fromResults(chunks[8], chunks[2], chunks[3], Double.parseDouble(chunks[4]), resultsType);
				
				/* Get the reigning match, if better, add */
				Match_fromResults otherMatch = bestResults.get(match.getSpectrumMD5());
				if (otherMatch == null) {
					bestResults.put(match.getSpectrumMD5(), match);
				} else {
					if (match.getScore() > otherMatch.getScore()) {
						bestResults.put(match.getSpectrumMD5(), match);
					}
				}
				
				/* at the end of it all, read in the new line */
				line = br.readLine();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	
	

}
