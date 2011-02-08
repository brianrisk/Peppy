package Reports;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import Peppy.Definitions;
import Peppy.Match;
import Peppy.MatchPTM;
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
		printHeader();
		Spectrum spectrum = match.getSpectrum();
		Peptide peptide = match.getPeptide();
		String acidString = peptide.getAcidSequenceString();
		printP("from protein: " + peptide.getProtein().getName());
		printP("spectum: " +match.getSpectrum().getFile().getName());
		printP("peptide: " + acidString);
		printP("peptide start: " + peptide.getStartIndex());
		printP("IMP value: " + match.getImpValue());
		
		//input spectrum image
		File spectrumFile = new File(destinationFile.getParentFile(), "spectrum.jpg");
		printP("<img src=\"spectrum.jpg\">");
		
		if (match.hasModification()) {
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
				printLI(acidString.charAt(i) + ": " + imp);
				if (imp < bestIMP) {
					bestIMP = imp;
					bestIndex = i;
					bestMatch = matchPTM;
				}
			}
			print("</ol>");
			
			//printing the spectrum match with the best mod location we just found
			try {
				SpectralVisualizerPTM.drawDeluxSpectrum(spectrum, peptide, spectrumFile, bestMatch.getDifference(), bestIndex);
			} catch (IOException e) {
				e.printStackTrace();
			}			
		} else {
			try {
				SpectralVisualizer.drawDeluxSpectrum(spectrum, peptide, spectrumFile);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		printFooter();
	}

}
