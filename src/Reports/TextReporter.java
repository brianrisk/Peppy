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
	
	
	/**
	 * @param matches
	 * @param spectra
	 * @param sequences
	 */
	public TextReporter(ArrayList<Match> matches,
			ArrayList<Spectrum> spectra, ArrayList<Sequence> sequences) {
		this.matches = matches;
		this.spectra = spectra;
		this.sequences = sequences;
	}


	public void generateFullReport() {
		File reportFile = new File(Properties.reportDirectory, "textReport.txt");
		try {
			//create our report directory
			Properties.reportDirectory.mkdirs();
			
			//set up our main index file
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(reportFile)));
			
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
			sb.append("Precursor");
			sb.append('\t');
			sb.append("E Value");
			sb.append('\t');
			sb.append("Peptide");
			sb.append('\t');
			if (Peppy.Properties.isSequenceFileDNA) {
				sb.append("Sequence File");
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
			
			pw.println(sb);
			
			
			//print rows
			for (Match match: matches) {;
				sb = new StringBuffer();
				sb.append(match.getSpectrum().getId());
				sb.append('\t');
				sb.append(match.getSpectrum().getMD5());
				sb.append('\t');
				sb.append(match.getSpectrum().getFile().getName());
				sb.append('\t');
				sb.append(match.getScore());
				sb.append('\t');
				sb.append(match.getSpectrum().getPrecursorMass());
				sb.append('\t');
				sb.append(match.getEValue());
				sb.append('\t');
				sb.append(match.getPeptide().getAcidSequence());
				sb.append('\t');
				if (Peppy.Properties.isSequenceFileDNA) {
					sb.append(match.getSequence().getSequenceFile().getName());
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
					sb.append(match.getPeptide().getProteinName());
				}
				sb.append('\t');
				sb.append(match.getRank());
				
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
