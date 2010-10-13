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
public class BEDReporter {
	
	ArrayList<Match> matches;
	ArrayList<Spectrum> spectra;
	ArrayList<Sequence> sequences;
	
	
	/**
	 * @param matches
	 * @param spectra
	 * @param sequences
	 */
	public BEDReporter(ArrayList<Match> matches,
			ArrayList<Spectrum> spectra, ArrayList<Sequence> sequences) {
		this.matches = matches;
		this.spectra = spectra;
		this.sequences = sequences;
	}


	public void generateFullReport() {
		File reportFile = new File(Properties.reportDirectory, "BEDReport.bed");
		try {
			//create our report directory
			Properties.reportDirectory.mkdirs();
			
			//set up our main index file
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(reportFile)));
			
			//sorting our matches by spectrum then score
			Match.setSortParameter(Match.SORT_BY_E_VALUE);
			Collections.sort(matches);
			
			//print the header
			pw.println("#chrom	chromStart	chromEnd	name	score");
			
			
			String chrom;
			String chromFileName;
			int chromStart;
			int chromEnd;
			String name;
			double rawScore;
			int score;
			for (Match match: matches) {
				//chrom
				chromFileName = match.getSequence().getSequenceFile().getName();
				chrom = chromFileName.substring(0, chromFileName.indexOf('.'));
				
				//chromStart
				chromStart = match.getPeptide().getStartIndex();
				
				//chromEnd
				if (match.getPeptide().isForward()) {
					chromEnd = chromStart + 3 * match.getPeptide().getAcidSequence().length();
				} else {
					chromEnd = chromStart - 3 * match.getPeptide().getAcidSequence().length();
				}
				
				//name
				name = match.getSpectrum().getFile().getName();
				
				//score
				rawScore = Math.log(match.getEValue());
				score = 0;
				if (rawScore < -0.1) {
					score = (int) rawScore * -20;
					if (score > 1000) score = 1000;
				}
				
				pw.println(chrom + "\t" + chromStart + "\t" + chromEnd + "\t" + name + "\t" + score);
				
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
