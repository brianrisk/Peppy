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
import Peppy.U;


/**
 * Okay, so you've got all of your results.  Now what?
 * I'll tell you now what.  You want to see them presented in
 * a nice, easy and intuitive manner.  That's what this
 * class does.
 * 
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class TextReporter {
	
	ArrayList<Match> matches;
	File reportDir;
	
	
	/**
	 * @param matches
	 * @param spectra
	 * @param sequence_DNAs
	 */
	public TextReporter(ArrayList<Match> matches, File reportDir) {
		this.matches = matches;
		this.reportDir = reportDir;
	}

	public void generateFullReport() {	
		reportDir.mkdirs();
		//set up our main index file
		File reportFile = new File(reportDir, "report.txt");
		try {	
			
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(reportFile)));
			
			/* CHANGE THIS WITH EACH ADJUSTMENT TO FILE FORMAT */
			pw.println("format version 19");
			
			if (Properties.isSequenceFileDNA) {
				pw.println("> analysis-type: nucleotide");
			} else {
				pw.println("> analysis-type: protein");
			}
			pw.println("> scoring method:" + Properties.scoringMethodName);
			
			//sorting our matches by spectrum then score
			Match.setSortParameter(Match.SORT_BY_SCORE);
			Collections.sort(matches);
			
			
			StringBuffer sb;
			//print header
			sb = new StringBuffer();
			sb.append("spectrumID");
			sb.append('\t');
			sb.append("fileLocus");
			sb.append('\t');
			sb.append("spectrumMD5");
			sb.append('\t');
			sb.append("FilePath");
			sb.append('\t');
			sb.append("score");
			sb.append('\t');
			sb.append("peptideMass");
			sb.append('\t');
			sb.append("PrecursorNeutralMass");
			sb.append('\t');
			sb.append("peptideSequence");
			sb.append('\t');
			sb.append("previousAminoAcid");
			sb.append('\t');
			sb.append("start");
			sb.append('\t');
			sb.append("stop");
			sb.append('\t');
			sb.append("SequenceName");
			if (Peppy.Properties.isSequenceFileDNA) {
				sb.append('\t');
				sb.append("INTRON-START");
				sb.append('\t');
				sb.append("INTRON-STOP");
				sb.append('\t');
				sb.append("Strand");
				sb.append('\t');
				sb.append("Is Spliced");
			}
			sb.append('\t');
			sb.append("RankCount");
			sb.append('\t');
			sb.append("Charge");	
			sb.append('\t');
			sb.append("inORF");
			sb.append('\t');
			sb.append("sizeOfORF");
			sb.append('\t');
			sb.append("Hydrophobic");
			sb.append('\t');
			sb.append("Hydrophilic");
			sb.append('\t');
			sb.append("isModified");
			sb.append('\t');
			sb.append("modMass");
			sb.append('\t');
			sb.append("modIndex");
			sb.append('\t');
			sb.append("modLocCertain");
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
