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
			pw.println("format version 9");
			
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
			sb.append("Spectrum ID");
			sb.append('\t');
			sb.append("MD5");
			sb.append('\t');
			sb.append("File Name");
			sb.append('\t');
			sb.append("Score");
			sb.append('\t');
			sb.append("Precursor M/Z");
			sb.append('\t');
			sb.append("Precursor Neutral Mass");
			sb.append('\t');
			sb.append("E Value");
			sb.append('\t');
			sb.append("Peptide Sequence");
			sb.append('\t');
			if (Peppy.Properties.isSequenceFileDNA) {
				sb.append("Sequence File");
				sb.append('\t');
				sb.append("Sequence Description");
				sb.append('\t');
				sb.append("START");
				sb.append('\t');
				sb.append("STOP");
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
			sb.append("Match Rank");
			sb.append('\t');
			sb.append("Rank Count");
			sb.append('\t');
			sb.append("Ion Count");	
			sb.append('\t');
			sb.append("Labeled");	
			sb.append('\t');
			sb.append("Charge");	
			sb.append('\t');
			sb.append("CleavageAcidCount");
			pw.println(sb);		
			
			//print rows
			for (Match match: matches) {
				sb = new StringBuffer();
				sb.append(match.getSpectrum().getId());
				sb.append('\t');
				sb.append(match.getSpectrum().getMD5());
				sb.append('\t');
				sb.append(match.getSpectrum().getFile().getName());
				sb.append('\t');
				sb.append(match.getScore());
				sb.append('\t');
				sb.append(match.getSpectrum().getPrecursorMZ());
				sb.append('\t');
				sb.append(match.getSpectrum().getMass());
				sb.append('\t');
				sb.append(match.getEValue());
				sb.append('\t');
				sb.append(match.getPeptide().getAcidSequenceString());
				sb.append('\t');
				if (Peppy.Properties.isSequenceFileDNA) {
					sb.append(match.getPeptide().getParentSequence().getSequenceFile().getName());
					sb.append('\t');
					sb.append(match.getPeptide().getProtein().getName());
					sb.append('\t');
					sb.append(match.getPeptide().getStartIndex());
					sb.append('\t');
					sb.append(match.getPeptide().getStopIndex());
					sb.append('\t');
					sb.append(match.getPeptide().getIntronStartIndex());
					sb.append('\t');
					sb.append(match.getPeptide().getIntronStopIndex());
					sb.append('\t');
					sb.append(match.getPeptide().isForward() ? "+" : "-");
					sb.append('\t');
					sb.append(match.getPeptide().isSpliced());
				} else {
					sb.append(match.getPeptide().getProtein().getName());
				}
				sb.append('\t');
				sb.append(match.rank);
				sb.append('\t');
				sb.append(match.repeatCount);
				sb.append('\t');
				sb.append(match.getIonMatchTally());
				sb.append('\t');
				sb.append(match.isIsotopeLabeled());
				sb.append('\t');
				sb.append(match.getSpectrum().getCharge());
				sb.append('\t');
				sb.append(match.getPeptide().getCleavageAcidCount());
				pw.println(sb);
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
