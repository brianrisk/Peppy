package Reports;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;

import Peppy.Match;
import Peppy.Properties;
import Peppy.Regions;
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
public class HTMLReporter {
	
	ArrayList<Match> matches;
	ArrayList<Spectrum> spectra;
	ArrayList<Sequence> sequences;
	File reportDir;
	
	
	/**
	 * @param matches
	 * @param spectra
	 * @param sequence_DNAs
	 */
	public HTMLReporter(ArrayList<Match> matches,
			ArrayList<Spectrum> spectra, ArrayList<Sequence> sequences, File reportDir) {
		this.matches = matches;
		this.spectra = spectra;
		this.sequences = sequences;
		this.reportDir = reportDir;
	}


	public void generateFullReport() {
		//create our report directory
		reportDir.mkdirs();
		File indexFile = new File(reportDir, "index" + Properties.reportWebSuffix);

		//sorting our matches by spectrum then score
		Match.setSortParameter(Match.SORT_BY_SPECTRUM_ID_THEN_SCORE);
		Collections.sort(matches);
		
		ArrayList<Match> bestMatches = new ArrayList<Match>();
		
		int spectrumID = -1; //an ID will never be -1
		int matchRank = 0;
		double scoreRatio;
		Match previousMatch = null;
		for (Match match: matches) {
			if (spectrumID != match.getSpectrum().getId()) {
				spectrumID = match.getSpectrum().getId();
				matchRank = 1;
				bestMatches.add(match);
			} else {
				matchRank++;
			}
			if (matchRank == 2) {
				scoreRatio = previousMatch.getScore() / match.getScore();
				previousMatch.setScoreRatio(scoreRatio);
			}
			previousMatch = match;
		}
		
		//sort our best matches
		Match.setSortParameter(Match.SORT_BY_IMP_VALUE);
		Collections.sort(bestMatches);
		
		NumberFormat nfDecimal = NumberFormat.getInstance();
		nfDecimal.setMaximumFractionDigits(4);
		NumberFormat nfPercent = NumberFormat.getPercentInstance();
		nfPercent.setMaximumFractionDigits(2);
		
		//limit how many we will display
		int maxDisplay = 10000;
		if (maxDisplay > bestMatches.size()) maxDisplay = bestMatches.size();
		
		HTMLPageMatches hpm = new HTMLPageMatches(bestMatches, maxDisplay, indexFile);
		hpm.makePage();

		/* create spectrum reports */
		if (Properties.generateSpectrumReport) {
			for (Match match: bestMatches) {
				generateSpectrumReport(match.getSpectrum());
			}
		}
		
		/* create region report */
		Regions regions = new Regions(matches, sequences, spectra);
		regions.createReport(reportDir);
		regions.clearRegions();

	}
	
	
	
	
	public void generateSpectrumReport(Spectrum spectrum) {
		//set up our files and folders
		File sequenceDirectory = new File(reportDir, "spectra");
		sequenceDirectory.mkdirs();
		File indexFile = new File(sequenceDirectory, spectrum.getId() + Properties.reportWebSuffix);
		
		//our report object
		HTMLPage page;
		

		page = new HTMLPageSpectrum(spectrum, matches, indexFile);

		page.makePage();
	}
	
	

	
	
}
