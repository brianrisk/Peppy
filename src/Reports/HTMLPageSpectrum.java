package Reports;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;

import Math.MassError;
import Peppy.Definitions;
import Peppy.Match;
import Peppy.Match_IMP_VariMod;
import Peppy.Matches;
import Peppy.Modification;
import Peppy.Peak;
import Peppy.Properties;
import Peppy.Spectrum;

public class HTMLPageSpectrum extends HTMLPage {

	ArrayList<Match> theseMatches;
	Spectrum spectrum;
	
	public HTMLPageSpectrum(Spectrum spectrum, ArrayList<Match> matches, File destinationFile) {
		super(destinationFile);
		this.spectrum = spectrum;
		theseMatches =  Matches.getMatchesWithSpectrum(spectrum, matches);
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
					spectrumScript.append(match_IMP_VariMod.getMoificationdMass());
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
		spectrumScript.append("<script src=\"http://peppyresearch.com/js/processing-1.1.0.js\"></script>");
		spectrumScript.append("<script src=\"http://peppyresearch.com/spectrumvisualizer/psv-control.js\"></script>");
		
		//print header
		printHeader("Spectrum report for " + spectrum.getFile().getName(), spectrumScript.toString());
		
		printP("<ul><li>mass: " + spectrum.getMass() + "<li>charge: " + spectrum.getCharge() + "<li>m/z: " + spectrum.getPrecursorMZ() + "</ul>");
		
		//spectrum
		printP("<canvas data-processing-sources=\"http://peppyresearch.com/spectrumvisualizer/PeppySpectrumVisualizer.pjs\" id=\"spectrum\" width=\"800\" height=\"310\"></canvas>");
//		printP("<canvas data-processing-sources=\"ionMatchVisualizer.pjs\" id=\"spectrum\" width=\"800\" height=\"310\"></canvas>");
		
		/* modifications */
		if (theseMatches.get(0).hasModification()) {
			printP("known modifications which are close in mass to the observed modification:");
			print("<ul>");
			Match_IMP_VariMod modifiedMatch = (Match_IMP_VariMod) theseMatches.get(0);
			double fragmentTolerance = MassError.getDaltonError(Properties.fragmentTolerance, 1000);
			for (Modification mod: Definitions.modifications) {
				if (Math.abs(modifiedMatch.getMoificationdMass() - mod.getMonoMass()) <= fragmentTolerance) {
					print("<li>" + mod.getDescription() + "</li>");
				}
			}
			print("</ul>");
		}
		
		
		//Our table
		print("<table class=\"sortable\" id=\"box-table-a\" width=\"95%\">");
		printTR();
		printTH("UCSC");
		printTH("Google");
		printTH("peptide");
		printTH("sequence");
		printTH("indicies");
		printTH("F");
		if (Properties.useSpliceVariants) {
			printTH("S");
		}
		printTH("score");
		printTH("ions");
//		if (Properties.searchModifications) {
			printTH("has mod");
			printTH("mod index");
			printTH("mod mass");
//		}
		/* print all the rows */
		for(Match match: theseMatches) {
			printTableRow(match);
		}
		print("</table>");
		
		//histogram
//		printH2("E value histogram for spectrum " + spectrum.getId());
//		File histogramFile = new File(destinationFile.getParent(), spectrum.getId() + "-hist.jpg");
//		try {
//			HistogramVisualizer.drawHistogram(spectrum.getEValueCalculator().getSmoothedHistogram(), 300, 300, histogramFile);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//		printP("<img src=\"" + histogramFile.getName() + "\">");

		printFooter();
	}
	
	private void printTableRow(Match match) {
		printTR();
		
		/* UCSC link */
		String link = UCSC.getLink(match);
		printTD("(<a href=\"" +link + "\">UCSC</a>)");
		
		/* google link */
		link = "http://www.google.com/search?hl=en&q=" + match.getPeptide().getAcidSequenceString();
		printTD("(<a href=\"" +link + "\">Google</a>)");
		
		
		//peptide sequence
		StringBuffer peptideLine = new StringBuffer();
		peptideLine.append("<a href=\"\" class=\"spectrumTrigger\" onClick=\"javascript:changePeptide('");
		peptideLine.append(match.getPeptide().getAcidSequenceString());
		peptideLine.append("');} return false;\">");
		peptideLine.append(match.getPeptide().getAcidSequenceString());
		peptideLine.append("</a>");
		printTD(peptideLine.toString());
		
		//sequence name
		if (Properties.useSpliceVariants) {
			printTD("NULL");
		} else {
			printTD(match.getPeptide().getProtein().getName());
		}
		//start / stop
		printTD( match.getPeptide().getStartIndex() + ", " + match.getPeptide().getStopIndex());
		
		//forward/reverse
		if (match.getPeptide().isForward()) {printTD("+");}
		else {printTD("-");}
		
		//isSpliced
		if (Properties.useSpliceVariants) {
			if (match.getPeptide().isSpliced()) {printTD("Y");}
			else {printTD("N");}
		}
		
		//score
		printTD("" + match.getScore());
		
		//ion matdch tally
		printTD("" + match.getIonMatchTally());
		
//		if (Properties.searchModifications) {
			printTD("" + match.hasModification());
			if (match.hasModification()) {
				printTD((match.getModificationIndex() + 1) + " (" + match.getPeptide().getAcidSequenceString().charAt(match.getModificationIndex()) + ")");
				printTD("" + match.getMoificationdMass());
			} else {
				printTD("");
				printTD("");
			}
			
			
//		}
	}
	

}
