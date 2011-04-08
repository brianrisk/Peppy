package Reports;

import java.io.File;
import java.util.ArrayList;

import Peppy.Definitions;
import Peppy.Match;
import Peppy.Match_IMP_VariMod;
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
		//TODO this only accounts for vari_mods when this could be a multi mod
		if (match.hasModification()) {
			Match_IMP_VariMod match_IMP_VariMod = (Match_IMP_VariMod) match;
			spectrumScript.append("var modifications = [");
			for (int i = 0; i < peptide.getLength(); i++) {
				if (i == match_IMP_VariMod.getModificationIndex()) {
					spectrumScript.append(match_IMP_VariMod.getModificationMass());
				} else {
					spectrumScript.append("0");
				}
				if (i <  peptide.getLength() - 1) spectrumScript.append(", ");
			}
			spectrumScript.append("];");
		} else {
			spectrumScript.append("var modifications = [];");
		}
		spectrumScript.append("var acidSequence = '" + acidString + "';");
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
		Modification [] modifications = match.getModifications();
		for (int i = 0; i < modifications.length; i++) {
			spectrumScript.append(modifications[i].getMonoMass());
			if (i < peaks.size() - 1) spectrumScript.append(", ");
		}
		spectrumScript.append("];");
		

		spectrumScript.append("function changePeptide() {");
		spectrumScript.append("	acidSequence = document.forms['acidSequence'].elements['sequence'].value.toUpperCase();");
		spectrumScript.append("	Processing.getInstanceById('spectrum').markMatchingIons();	");
		spectrumScript.append("}");
		spectrumScript.append("");
		spectrumScript.append("var peakDifferenceThreshold = 0.5;");
		spectrumScript.append("function changePeakDifferenceThreshold() {");
		spectrumScript.append("	peakDifferenceThreshold = parseFloat(document.forms['peakThreshold'].elements['value'].value);");
		spectrumScript.append("	Processing.getInstanceById('spectrum').markMatchingIons();");
		spectrumScript.append("}");
		spectrumScript.append("");
		spectrumScript.append("//alternates between scaled y and not scaled");
		spectrumScript.append("var scaleYAxis = false;");
		spectrumScript.append("function toggleYAxisDisplay() {");
		spectrumScript.append("	if (!scaleYAxis) {");
		spectrumScript.append("		document.forms['flipScale'].elements['submit'].value = 'un-scale y-axis';");
		spectrumScript.append("	} else {");
		spectrumScript.append("		document.forms['flipScale'].elements['submit'].value = 'scale y-axis';");
		spectrumScript.append("	}");
		spectrumScript.append("	scaleYAxis = !scaleYAxis;");
		spectrumScript.append("}");
		spectrumScript.append("	");
		spectrumScript.append("");
		spectrumScript.append("//changes if we see masses for each peak");
		spectrumScript.append("var displayMasses = false;");
		spectrumScript.append("function toggleMassDisplay() {");
		spectrumScript.append("	if (!displayMasses) {");
		spectrumScript.append("		document.forms['flipMasses'].elements['submit'].value = 'hide masses';");
		spectrumScript.append("	} else {");
		spectrumScript.append("		document.forms['flipMasses'].elements['submit'].value = 'show masses';");
		spectrumScript.append("	}");
		spectrumScript.append("	displayMasses = !displayMasses;");
		spectrumScript.append("}");
		
		spectrumScript.append("</script>");
		spectrumScript.append("<script src=\"http://peppyresearch.com/js/processing-1.0.0.js\"></script>");
//		spectrumScript.append("<script src=\"../../processing-1.0.0.js\"></script>");
		
		//print header
		printHeader("Spectrum report for " + spectrum.getFile().getName(), spectrumScript.toString());
		
		//spectrum
		printP("<canvas data-processing-sources=\"http://peppyresearch.com/spectrumvisualizer/PeppySpectrumVisualizer.pjs\" id=\"spectrum\" width=\"800\" height=\"310\"></canvas>");
//		printP("<canvas data-processing-sources=\"../../ionMatchVisualizer.pjs\" id=\"spectrum\" width=\"800\" height=\"310\"></canvas>");
		
//		printP("from protein: " + peptide.getProtein().getName());
//		printP("peptide: " + acidString);
//		printP("peptide start: " + peptide.getStartIndex());
//		printP("IMP value: " + match.getImpValue());
		
		if (match.hasModification()) {
			//this is a Modification match, make it so
			Match_IMP_VariMod match_IMP_VariMod = (Match_IMP_VariMod) match;
			
			//details table
			print("<table valign=\"top\" width=\"95%\">");
			printTR();
			print("<td>");
			
			printH2("Modification properties");
			printP("mass difference: " + match_IMP_VariMod.getModificationMass());
			
			//print the probable modifications
			ArrayList<Modification> potentialModifications = new ArrayList<Modification>(); 
			for (Modification pm: Definitions.modifications) {
				if (Math.abs(match_IMP_VariMod.getModificationMass() - pm.getMonoMass()) < 0.1) {
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
				Match_IMP_VariMod indexModificationMatch = new Match_IMP_VariMod(spectrum, peptide);
				double indexIMP = indexModificationMatch.calculateIMP(match_IMP_VariMod.getModificationMass(), i);
				
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
