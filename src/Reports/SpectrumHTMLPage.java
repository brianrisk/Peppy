package Reports;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;

import Peppy.Match;
import Peppy.Match_IMP_VariMod;
import Peppy.Peak;
import Peppy.Spectrum;

public class SpectrumHTMLPage extends HTMLPage {

	ArrayList<Match> theseMatches;
	Spectrum spectrum;
	
	public SpectrumHTMLPage(Spectrum spectrum, ArrayList<Match> matches, File destinationFile) {
		super(destinationFile);
		this.spectrum = spectrum;
		theseMatches =  MatchSearches.getMatchesWithSpectrum(spectrum, matches);
		Match.setSortParameter(Match.SORT_BY_SCORE);
		Collections.sort(theseMatches);
	}

	@Override
	public void makePage() {
		//build spectrum script for header
		StringBuffer spectrumScript = new StringBuffer();
		spectrumScript.append("<script type=\"text/javascript\">");
		spectrumScript.append("var acidSequence = '");
		spectrumScript.append(theseMatches.get(0).getPeptide().getAcidSequenceString());
		spectrumScript.append("';");
		spectrumScript.append("var spectrumMass = " + spectrum.getMass() + ";");
		spectrumScript.append("var spectrumMaxIntensity = " + spectrum.getMaxIntensity() + ";");
		spectrumScript.append("var peakMasses = [");
		ArrayList<Peak> peaks = spectrum.getPeaks();
		for (int i = 0; i < peaks.size(); i++) {
			Peak peak = peaks.get(i);
			spectrumScript.append(peak.getMass());
			if (i < peaks.size() - 1) spectrumScript.append(", ");
		}
		spectrumScript.append("];");
		spectrumScript.append("var peakIntensities = [");
		for (int i = 0; i < peaks.size(); i++) {
			Peak peak = peaks.get(i);
			spectrumScript.append(peak.getIntensity());
			if (i < peaks.size() - 1) spectrumScript.append(", ");
		}
		spectrumScript.append("];");
		
		if (theseMatches.get(0).hasModification()) {
			Match_IMP_VariMod match_IMP_VariMod = (Match_IMP_VariMod) theseMatches.get(0);
			spectrumScript.append("var modifications = [");
			for (int i = 0; i < theseMatches.get(0).getPeptide().getLength(); i++) {
				if (i == match_IMP_VariMod.getModificationIndex()) {
					spectrumScript.append(match_IMP_VariMod.getModificationMass());
				} else {
					spectrumScript.append("0");
				}
				if (i <  theseMatches.get(0).getPeptide().getLength() - 1) spectrumScript.append(", ");
			}
			spectrumScript.append("];");
		} else {
			spectrumScript.append("var modifications = [];");
		}
		

		
		spectrumScript.append("</script>");
		spectrumScript.append("<script src=\"http://peppyresearch.com/js/processing-1.0.0.js\"></script>");
		spectrumScript.append("<script src=\"http://peppyresearch.com/spectrumvisualizer/pvs-control.js\"></script>");
//		spectrumScript.append("<script src=\"processing-1.0.0.js\"></script>");
		
		//print header
		printHeader("Spectrum report for " + spectrum.getFile().getName(), spectrumScript.toString());
		
		//histogram
//		printH2("E value histogram for spectrum " + spectrum.getId());
//		File histogramFile = new File(destinationFile.getParent(), spectrum.getId() + "-hist.jpg");
//		try {
//			HistogramVisualizer.drawHistogram(spectrum.getEValueCalculator().getSmoothedHistogram(), 300, 300, histogramFile);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//		printP("<img src=\"" + histogramFile.getName() + "\">");
		
		//spectrum
		printP("<canvas data-processing-sources=\"http://peppyresearch.com/spectrumvisualizer/PeppySpectrumVisualizer.pjs\" id=\"spectrum\" width=\"800\" height=\"310\"></canvas>");
//		printP("<canvas data-processing-sources=\"ionMatchVisualizer.pjs\" id=\"spectrum\" width=\"800\" height=\"310\"></canvas>");
		
		//Our table
		print("<table class=\"sortable\" id=\"box-table-a\" width=\"95%\">");
		printTR();
		printTH("peptide");
		printTH("database");
		printTH("indicies");
		printTH("F");
		printTH("S");
		printTH("score");
		printTH("ions");
		printTH("E value");
		for(Match match: theseMatches) {
			printTableRow(match);
		}
		print("</table>");

		printFooter();
	}
	
	private void printTableRow(Match match) {
		printTR();
		
		//peptide sequence
		StringBuffer peptideLine = new StringBuffer();
		peptideLine.append("<a href=\"\" class=\"spectrumTrigger\" onClick=\"javascript:changePeptide('");
		peptideLine.append(match.getPeptide().getAcidSequenceString());
		peptideLine.append("');} return false;\">");
		peptideLine.append(match.getPeptide().getAcidSequenceString());
		peptideLine.append("</a>");
		printTD(peptideLine.toString());
		
		//sequence name
		if (match.getPeptide().getParentSequence() != null) {
			printTD(match.getPeptide().getParentSequence().getSequenceFile().getName());
		}
		
		//start / stop
		printTD( match.getPeptide().getStartIndex() + ", " + match.getPeptide().getStopIndex());
		
		//forward/reverse
		if (match.getPeptide().isForward()) {printTD("+");}
		else {printTD("-");}
		
		//isSpliced
		if (match.getPeptide().isSpliced()) {printTD("Y");}
		else {printTD("N");}
		
		//score
		printTD("" + match.getScore());
		
		//ion matdch tally
		printTD("" + match.getIonMatchTally());
		
		//E value
		printTD("" + match.getEValue());
	}
	

}
