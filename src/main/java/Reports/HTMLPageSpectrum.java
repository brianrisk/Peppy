package Reports;

import Math.MassError;
import Peppy.*;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
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
		spectrumScript.append("var spectrumMass = " + spectrum.getMass()  + ";");
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
		spectrumScript.append("<script src=\"http://geneffects.com/files/peppy/js/processing-1.1.0.js\"></script>");
//		spectrumScript.append("<script src=\"http://rrcs-98-101-157-178.midsouth.biz.rr.com/~risk2/js/processing.js\"></script>");
		spectrumScript.append("<script src=\"http://geneffects.com/files/peppy/spectrumvisualizer/psv-control.js\"></script>");
//		spectrumScript.append("<script src=\"http://rrcs-98-101-157-178.midsouth.biz.rr.com/~risk2/spectrumvisualizer/psv-control.js\"></script>");

		
		//print header
		printHeader("Spectrum report for " + spectrum.getFile().getName(), spectrumScript.toString());
		
		printP("<ul><li>neutral mass: " + spectrum.getMass() + "<li>charge: " + spectrum.getCharge() + "</ul>");
		
		//spectrum
		printP("<canvas data-processing-sources=\"http://geneffects.com/files/peppy/spectrumvisualizer/PeppySpectrumVisualizer.pjs\" id=\"spectrum\" width=\"800\" height=\"310\"></canvas>");
//		printP("<canvas data-processing-sources=\"http://rrcs-98-101-157-178.midsouth.biz.rr.com/~risk2/spectrumvisualizer/PeppySpectrumVisualizer.pjs\" id=\"spectrum\" width=\"800\" height=\"310\"></canvas>");
		
//		printP("<canvas data-processing-sources=\"ionMatchVisualizer.pjs\" id=\"spectrum\" width=\"800\" height=\"310\"></canvas>");
		
		/* modifications */
		if (theseMatches.get(0).hasModification()) {
			printP("known modifications which are close in mass to the observed modification:");
			print("<ul>");
			Match_IMP_VariMod modifiedMatch = (Match_IMP_VariMod) theseMatches.get(0);
			double fragmentTolerance = MassError.getDaltonError(Properties.fragmentTolerance, spectrum.getMass());
			for (ModificationEntry mod: Definitions.modificationEntries) {
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
		printTH("score");
		printTH("has mod");
		printTH("mod index");
		printTH("mod mass");
		
		/* print all the rows */
		for(Match match: theseMatches) {
			if (!(match instanceof Match_Blank)) {
				printTableRow(match);
			}
		}
		print("</table>");
		

		printFooter();
	}
	
	private void printTableRow(Match match) {
		printTR();
		
		/* UCSC link */
		String link = "";
		link = UCSC.getLink(match);
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
		printTD(match.getPeptide().getProtein().getName());
		//start / stop
		printTD( match.getPeptide().getStartIndex() + ", " + match.getPeptide().getStopIndex());
		
		//forward/reverse
		if (match.getPeptide().isForward()) {printTD("+");}
		else {printTD("-");}
		
		//score
		printTD("" + match.getScore());
		

		printTD("" + match.hasModification());
		if (match.hasModification()) {
			printTD((match.getModificationIndex() + 1) + " (" + match.getPeptide().getAcidSequenceString().charAt(match.getModificationIndex()) + ")");
			printTD("" + match.getMoificationdMass());
		} else {
			printTD("");
			printTD("");
		}

			

	}
	

}
