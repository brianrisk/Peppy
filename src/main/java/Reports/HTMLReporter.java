package Reports;

import Peppy.Match;
import Peppy.Properties;
import Peppy.Spectrum;

import java.io.File;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;


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
public class HTMLReporter {
	
	ArrayList<Match> matches;
	File reportDir;
	
	
	/**
	 * @param matches
	 * @param spectra
	 * @param sequence_DNAs
	 */
	public HTMLReporter(ArrayList<Match> matches, File reportDir) {
		this.matches = matches;
		this.reportDir = reportDir;
	}


	public void generateFullReport() {
		//create our report directory
		reportDir.mkdirs();
		File indexFile = new File(reportDir, "index" + Properties.reportWebSuffix);

		/* create a list of best matches */
		ArrayList<Match> bestMatches = new ArrayList<Match>();
		Match.setSortParameter(Match.SORT_BY_SPECTRUM_ID_THEN_SCORE);
		Collections.sort(matches);
		int spectrumID = -1; //an ID will never be -1
		for (Match match: matches) {
			if (spectrumID != match.getSpectrum().getId()) {
				spectrumID = match.getSpectrum().getId();
				bestMatches.add(match);
			}
		}
		
		//sort our best matches
		Match.setSortParameter(Match.SORT_BY_SCORE);
		Collections.sort(bestMatches);
		
		NumberFormat nfDecimal = NumberFormat.getInstance();
		nfDecimal.setMaximumFractionDigits(4);
		NumberFormat nfPercent = NumberFormat.getPercentInstance();
		nfPercent.setMaximumFractionDigits(2);
		
		//limit how many we will display
		int maxDisplay = 1000;
		if (maxDisplay > bestMatches.size()) maxDisplay = bestMatches.size();
		
		HTMLPageMatches hpm = new HTMLPageMatches(bestMatches, maxDisplay, indexFile);
		hpm.makePage();

		/* create spectrum reports */
		for (Match match: bestMatches) {
			generateSpectrumReport(match.getSpectrum());
		}
		

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
