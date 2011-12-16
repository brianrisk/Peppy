package Reports;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;

import Peppy.Match;
import Peppy.Properties;
import Peppy.Sequence;
import Peppy.Spectrum;
import Utilities.U;


/**
 * Okay, so you've got all of your results.  Now what?
 * I'll tell you now what.  You want to see them presented in
 * a nice, easy and intuitive manner.  That's what this
 * class does.
 * @author Brian Risk
 *
 */
public class TextReporter {
	
	ArrayList<Match> matches;
	ArrayList<Spectrum> spectra;
	ArrayList<Sequence> sequences;
	File reportDir;
	
	
	/**
	 * @param matches
	 * @param spectra
	 * @param sequence_DNAs
	 */
	public TextReporter(ArrayList<Match> matches,
			ArrayList<Spectrum> spectra, ArrayList<Sequence> sequences, File reportDir) {
		this.matches = matches;
		this.spectra = spectra;
		this.sequences = sequences;
		this.reportDir = reportDir;
	}

	public void generateFullReport() {	
		reportDir.mkdirs();
		//set up our main index file
		File reportFile = new File(reportDir, reportDir.getName() + "_report.txt");
		try {	
			
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(reportFile)));
			
			//CHANGE THIS WITH EACH ADJUSTMENT TO FILE FORMAT
			pw.println("format version 16");
			
			if (Properties.isSequenceFileDNA) {
				pw.println("> analysis-type: nucleotide");
			} else {
				pw.println("> analysis-type: protein");
			}
			pw.println("> scoring method:" + Properties.scoringMethodName);
			
			//sorting our matches by spectrum then score
			Match.setSortParameter(Match.SORT_BY_E_VALUE);
			Collections.sort(matches);
			
			
			StringBuffer sb;
			//print header
			sb = new StringBuffer();
			sb.append("ID");
			sb.append('\t');
			sb.append("SpectrumID");
			sb.append('\t');
			sb.append("MD5");
			sb.append('\t');
			sb.append("FileName");
			sb.append('\t');
			sb.append("Score");
			sb.append('\t');
			sb.append("PrecursorM/Z");
			sb.append('\t');
			sb.append("PrecursorNeutralMass");
			sb.append('\t');
			sb.append("E Value");
			sb.append('\t');
			sb.append("PeptideSequence");
			sb.append('\t');
			sb.append("START");
			sb.append('\t');
			sb.append("STOP");
			sb.append('\t');
			if (Peppy.Properties.isSequenceFileDNA) {
				sb.append("Sequence File");
				sb.append('\t');
				sb.append("Sequence Description");
				sb.append('\t');
				sb.append("INTRON-START");
				sb.append('\t');
				sb.append("INTRON-STOP");
				sb.append('\t');
				sb.append("Strand");
				sb.append('\t');
				sb.append("Is Spliced");
			} else {
				sb.append("Protein Name");
			}
			sb.append('\t');
			sb.append("MatchRank");
			sb.append('\t');
			sb.append("RankCount");
			sb.append('\t');
			sb.append("IonCount");	
			sb.append('\t');
			sb.append("Labeled");	
			sb.append('\t');
			sb.append("Charge");	
			sb.append('\t');
			sb.append("CleavageAcidCount");
			sb.append('\t');
			sb.append("inORF");
			sb.append('\t');
			sb.append("Hydrophobic");
			sb.append('\t');
			sb.append("Hydrophilic");
			if (Properties.searchModifications) {
				sb.append('\t');
				sb.append("isModified");
				sb.append('\t');
				sb.append("modMass");
				sb.append('\t');
				sb.append("modIndex");
			}
			pw.println(sb);		
			
			//print rows
			for (Match match: matches) {
				pw.println(match.toString());
			}
			

			pw.flush();
			pw.close();

			
		} catch (FileNotFoundException e) {
			U.p("could not find file: " + reportFile.getName());
			e.printStackTrace();
		} catch (IOException e) {
			U.p("could not read file: " + reportFile.getName());
			e.printStackTrace();
		}
		
	}
	
	

}
