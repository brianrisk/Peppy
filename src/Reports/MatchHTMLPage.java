package Reports;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import Peppy.Definitions;
import Peppy.Match;
import Peppy.MatchPTM;
import Peppy.Peak;
import Peppy.Peptide;
import Peppy.ProteinModification;
import Peppy.Spectrum;
import SpectralVisualizer.SpectralVisualizer;
import SpectralVisualizer.SpectralVisualizerPTM;

public class MatchHTMLPage extends HTMLPage {
	
	Match match;
	public MatchHTMLPage(Match match, File destinationFile) {
		super(destinationFile);
		this.match = match;
	}
	
	@Override
	public void makePage() {
		
		//get useful variables
		Spectrum spectrum = match.getSpectrum();
		Peptide peptide = match.getPeptide();
		String acidString = peptide.getAcidSequenceString();
		
		//build spectrum script for header
		StringBuffer spectrumScript = new StringBuffer();
		spectrumScript.append("<script type=\"text/javascript\">");
		if (match.hasModification()) {
			spectrumScript.append("modificationMass = ");
			spectrumScript.append(match.getSpectrum().getMass() - match.getPeptide().getMass());
			spectrumScript.append(";");
		}
		spectrumScript.append(" var acidSequence = '");
		spectrumScript.append(acidString);
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
		spectrumScript.append("</script>");
		spectrumScript.append("<script src=\"processing-1.0.0.js\"></script>");
		
		//print header
		printHeader("Spectrum report for " + spectrum.getFile().getName(), spectrumScript.toString());
		
		//spectrum
		printP("<canvas data-processing-sources=\"ionMatchVisualizer.pjs\" id=\"spectrum\" width=\"800\" height=\"310\"></canvas>");
		
//		printP("from protein: " + peptide.getProtein().getName());
//		printP("peptide: " + acidString);
//		printP("peptide start: " + peptide.getStartIndex());
//		printP("IMP value: " + match.getImpValue());
		
		if (match.hasModification()) {
			//details table
			print("<table valign=\"top\" width=\"95%\">");
			printTR();
			print("<td>");
			
			MatchPTM matchPTM = (MatchPTM) match;
			printH2("Modification properties");
			printP("mass difference: " + matchPTM.getDifference());
			
			//print the probable modifications
			ArrayList<ProteinModification> potentialModifications = new ArrayList<ProteinModification>(); 
			ArrayList<ProteinModification> proteinModifications = Definitions.proteinModifications;
			for (ProteinModification pm: proteinModifications) {
				if (Math.abs(matchPTM.getDifference() - pm.getMonoMass()) < 0.1) {
					potentialModifications.add(pm);
				}
			}
			if (potentialModifications.size() > 0) {
				printH3("The modification might be:");
				print("<table>");
				for (ProteinModification pm: potentialModifications) {
					print("<tr>");
					printTD(pm.getDescription());
					printTD("" + pm.getMonoMass());
				}
				print("</table>");
			} else {
				printP("no known single modification is close to the observed mass difference.");
			}
			
			
			print("</td>");
			print("<td>");
			
			//finding and printing best modification locations
			double imp;
			double bestIMP = Double.MAX_VALUE;
			MatchPTM bestMatch = null;
			int bestIndex = 0;
			printH3("Modification location scores:");
			print("<ol>");
			for (int i= 0; i < acidString.length(); i++) {
				matchPTM = new MatchPTM(spectrum, peptide);
				imp = matchPTM.calculateIMP(matchPTM.getDifference(), i);
				
				//building the link
				StringBuffer peptideLine = new StringBuffer();
				peptideLine.append("<a href=\"\" onMouseOver=\"javascript:modificationIndex='");
				peptideLine.append(i);
				peptideLine.append("'; return false;\">");
				
				//print out the acid string, with bold modification
				for (int j = 0; j <i; j++) {
					peptideLine.append(acidString.charAt(j));
				}
				peptideLine.append("<b>");
				peptideLine.append(acidString.charAt(i));
				peptideLine.append("</b>");
				for (int j = i+1; j <acidString.length(); j++) {
					peptideLine.append(acidString.charAt(j));
				}
				
				peptideLine.append(": " + imp);
				peptideLine.append("</a>");
				printLI(peptideLine.toString());
				if (imp < bestIMP) {
					bestIMP = imp;
					bestIndex = i;
					bestMatch = matchPTM;
				}
			}
			print("</ol>");
			

			print("</td>");
			print("</table>");
		}
		
		printFooter();
	}

}
