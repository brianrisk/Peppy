package Reports;

import java.io.File;
import java.util.ArrayList;

import Peppy.Definitions;
import Peppy.Match;
import Peppy.MatchPTM;
import Peppy.Modification;
import Peppy.Peak;
import Peppy.Peptide;
import Peppy.Spectrum;

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
			MatchPTM matchPTM = (MatchPTM) match;
			//modification index
			spectrumScript.append("modificationIndex = ");
			spectrumScript.append(matchPTM.getModificationIndex());
			spectrumScript.append(";");
			
			//modification index
			spectrumScript.append("modificationMass = ");
			spectrumScript.append(matchPTM.getDifference());
			spectrumScript.append(";");
		} else {
			spectrumScript.append("modificationIndex = -1; modificationMass = 0;");
		}
		spectrumScript.append("var lock = false;");
		spectrumScript.append("var acidSequence = '");
		spectrumScript.append(acidString);
		spectrumScript.append("';");
		spectrumScript.append("var spectrumMass = " + spectrum.getMass() + ";");
		spectrumScript.append("var spectrumMaxIntensity = " + spectrum.getMaxIntensity() + ";");
		
		//printing out spectrum's peak masses
		spectrumScript.append("var peakMasses = [");
		ArrayList<Peak> peaks = spectrum.getPeaks();
		for (int i = 0; i < peaks.size(); i++) {
			Peak peak = peaks.get(i);
			spectrumScript.append(peak.getMass());
			if (i < peaks.size() - 1) spectrumScript.append(", ");
		}
		spectrumScript.append("];");
		
		//printing our spectrum's peak intensities
		spectrumScript.append("var peakIntensities = [");
		for (int i = 0; i < peaks.size(); i++) {
			Peak peak = peaks.get(i);
			spectrumScript.append(peak.getIntensity());
			if (i < peaks.size() - 1) spectrumScript.append(", ");
		}
		spectrumScript.append("];");
		
		//print the modifications array
		spectrumScript.append("];");
		spectrumScript.append("var peakIntensities = [");
		ArrayList<Modification> modifications = match.getModifications();
		for (int i = 0; i < modifications.size(); i++) {
			Modification modification = modifications.get(i);
			spectrumScript.append(modification.getMonoMass());
			if (i < peaks.size() - 1) spectrumScript.append(", ");
		}
		spectrumScript.append("];");
		
		spectrumScript.append("</script>");
		spectrumScript.append("<script src=\"http://proteomics.me/resources/processing-1.0.0.js\"></script>");
//		spectrumScript.append("<script src=\"../../processing-1.0.0.js\"></script>");
		
		//print header
		printHeader("Spectrum report for " + spectrum.getFile().getName(), spectrumScript.toString());
		
		//spectrum
		printP("<canvas data-processing-sources=\"http://proteomics.me/resources/ionMatchVisualizer.pjs\" id=\"spectrum\" width=\"800\" height=\"310\"></canvas>");
//		printP("<canvas data-processing-sources=\"../../ionMatchVisualizer.pjs\" id=\"spectrum\" width=\"800\" height=\"310\"></canvas>");
		
//		printP("from protein: " + peptide.getProtein().getName());
//		printP("peptide: " + acidString);
//		printP("peptide start: " + peptide.getStartIndex());
//		printP("IMP value: " + match.getImpValue());
		
		if (match.hasModification()) {
			//this is a Modification match, make it so
			MatchPTM matchPTM = (MatchPTM) match;
			
			//details table
			print("<table valign=\"top\" width=\"95%\">");
			printTR();
			print("<td>");
			
			printH2("Modification properties");
			printP("mass difference: " + matchPTM.getDifference());
			
			//print the probable modifications
			ArrayList<Modification> potentialModifications = new ArrayList<Modification>(); 
			ArrayList<Modification> modifications = Definitions.modifications;
			for (Modification pm: modifications) {
				if (Math.abs(matchPTM.getDifference() - pm.getMonoMass()) < 0.1) {
					potentialModifications.add(pm);
				}
			}
			if (potentialModifications.size() > 0) {
				printH3("The modification might be:");
				print("<table>");
				for (Modification pm: potentialModifications) {
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
			
			//printing the list of modification points
			printH3("Modification location scores:");
			print("<ol>");
			//printing the modification possibility links
			
			for (int i= 0; i < acidString.length(); i++) {
				MatchPTM indexModificationMatch = new MatchPTM(spectrum, peptide);
				double indexIMP = indexModificationMatch.calculateIMP(matchPTM.getDifference(), i);
				
				//building the link
				StringBuffer modificationLink = new StringBuffer();
				if (indexIMP != match.getImpValue()) {
					modificationLink.append("<a href=\"\" class=\"spectrumTrigger\" onClick=\"javascript:lock=!lock;return false;\" onMouseOver=\"javascript:if (!lock) {modificationIndex='");
				} else {
					modificationLink.append("<a href=\"\" class=\"spectrumTrigger bestTrigger\" onClick=\"javascript:lock=!lock;return false;\" onMouseOver=\"if (!lock) {javascript:modificationIndex='");
				}
				modificationLink.append(i);
				modificationLink.append("';} return false;\">");
				
				//print out the acid string, with bold modification
				for (int j = 0; j <i; j++) {
					modificationLink.append(acidString.charAt(j));
				}
				modificationLink.append("<b>");
				modificationLink.append(acidString.charAt(i));
				modificationLink.append("</b>");
				for (int j = i+1; j <acidString.length(); j++) {
					modificationLink.append(acidString.charAt(j));
				}
				
				modificationLink.append(": " + indexIMP);
				modificationLink.append("</a>");
				printLI(modificationLink.toString());
			}
			print("</ol>");
			

			print("</td>");
			print("</table>");
		}
		
		printFooter();
	}

}
