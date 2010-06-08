package Reports;

import java.io.*;
import java.util.*;
import Peppy.*;
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
	
	ArrayList<SpectrumPeptideMatch> matches;
	ArrayList<Spectrum> spectra;
	ArrayList<Sequence> sequences;
	
	
	/**
	 * @param matches
	 * @param spectra
	 * @param sequences
	 */
	public TextReporter(ArrayList<SpectrumPeptideMatch> matches,
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
			SpectrumPeptideMatch.setSortParameter(SpectrumPeptideMatch.SORT_BY_E_VALUE);
			Collections.sort(matches);
			
			
			StringBuffer sb;
			for (SpectrumPeptideMatch match: matches) {
//				if (match.getScoreHMM() < 1) {
					
					sb = new StringBuffer();
					sb.append(match.getSpectrum().getId());
//					sb.append('\t');
//					sb.append(match.getSequence().getSequenceFile().getName());
	//				sb.append('\t');
	//				sb.append(match.getPeptide().getIndex());
//					sb.append('\t');
//					sb.append(match.getScoreHMM());
					sb.append('\t');
					sb.append(match.getScoreMSMSFit());
//					sb.append('\t');
//					sb.append(match.getMSMSFitScoreRatio());
//					sb.append('\t');
//					sb.append(match.getMSMSFitRank());
//					sb.append('\t');
//					sb.append(match.getSpectrum().getPrecursorMass());
//					sb.append('\t');
//					sb.append(match.getPeptide().getMass());
//					sb.append('\t');
//					sb.append(Math.abs(match.getPeptide().getMass() - match.getSpectrum().getPrecursorMass()));
					sb.append('\t');
					sb.append(match.getEValue());
					sb.append('\t');
					sb.append(match.getPeptide().getAcidSequence());
					sb.append('\t');
					sb.append(match.getSpectrum().getFile().getName());
					pw.println(sb);
//				}
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
