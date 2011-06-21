package Reports;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

import Peppy.Match;
import Peppy.Match_IMP_MultiMod;
import Peppy.Modification;
import Peppy.Peak;
import Peppy.Spectrum;

public class SpectrumMultiMod extends HTMLPage {

	ArrayList<Match> theseMatches;
	Spectrum spectrum;
	Match_IMP_MultiMod bestMatch;
	
	public SpectrumMultiMod(Spectrum spectrum, ArrayList<Match> matches, File destinationFile) {
		super(destinationFile);
		this.spectrum = spectrum;
		theseMatches =  MatchSearches.getMatchesWithSpectrum(spectrum, matches);
		Match.setSortParameter(Match.SORT_BY_SCORE);
		Collections.sort(theseMatches);
		bestMatch = (Match_IMP_MultiMod) theseMatches.get(0);
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
		
		//printing the peak mass array
		spectrumScript.append("var peakMasses = [");
		ArrayList<Peak> peaks = spectrum.getPeaks();
		for (int i = 0; i < peaks.size(); i++) {
			Peak peak = peaks.get(i);
			spectrumScript.append(peak.getMass());
			if (i < peaks.size() - 1) spectrumScript.append(", ");
		}
		spectrumScript.append("];");
		
		//printing the peak intensity array
		spectrumScript.append("var peakIntensities = [");
		for (int i = 0; i < peaks.size(); i++) {
			Peak peak = peaks.get(i);
			spectrumScript.append(peak.getIntensity());
			if (i < peaks.size() - 1) spectrumScript.append(", ");
		}
		spectrumScript.append("];");
		
		//printing the modifications array
		Modification [] bestModifications = bestMatch.getModifications();
		spectrumScript.append("var modifications = [");
		for (int i = 0; i < bestModifications.length; i++) {
			spectrumScript.append(bestModifications[i].getMonoMass());
			if (i < bestModifications.length - 1) spectrumScript.append(", ");
		}
		spectrumScript.append("];");
		
		//ending the script
		spectrumScript.append("</script>");
		spectrumScript.append("<script src=\"http://proteomics.me/resources/processing-1.1.0.js\"></script>");
		
		//print header
		printHeader("Spectrum report for " + spectrum.getFile().getName(), spectrumScript.toString());
		
		//histogram
		printH2("E value histogram for spectrum " + spectrum.getId());
		File histogramFile = new File(destinationFile.getParent(), spectrum.getId() + "-hist.jpg");
		try {
			HistogramVisualizer.drawHistogram(spectrum.getEValueCalculator().getSmoothedHistogram(), 300, 300, histogramFile);
		} catch (IOException e) {
			e.printStackTrace();
		}
		printP("<img src=\"" + histogramFile.getName() + "\">");
		
		//spectrum
		printP("<canvas data-processing-sources=\"http://proteomics.me/resources/ionMatchVisualizerMultiMod.pjs\" id=\"spectrum\" width=\"800\" height=\"310\"></canvas>");
//		printP("<canvas data-processing-sources=\"ionMatchVisualizer.pjs\" id=\"spectrum\" width=\"800\" height=\"310\"></canvas>");
		
		//a table which shows the mods of the best match
		print("<table>");
		String acid = bestMatch.getPeptide().getAcidSequenceString();
		for (int i = 0; i < bestModifications.length; i++) {
			printTR();
			printTD("" + acid.charAt(i));
			printTD(bestModifications[i].getDescription());
		}
		print("</table>");
		
		//peptide sequence
		printP("acid: " + acid);
		
		//sequence name
//		printP("sequence: " + bestMatch.getPeptide().getParentSequence().getSequenceFile().getName());
		
		//start / stop
		printP("indicies: " +  bestMatch.getPeptide().getStartIndex() + ", " + bestMatch.getPeptide().getStopIndex());
		
		//forward/reverse
		printP("is forward: " + bestMatch.getPeptide().isForward());
		
		//isSpliced
		printP("is spliced: " + bestMatch.getPeptide().isSpliced());
		
		//score
		printP("score: " + bestMatch.getScore());
		
		//ion match tally
		printP("ion match tally: " + bestMatch.getIonMatchTally());
		
		//E value
		printP("e value: " + bestMatch.getEValue());
		
		//a table of all matches
		printH3("All matches to this spectrum");
		print("<table>");
		for (Match matchGeneric: theseMatches) {
			Match_IMP_MultiMod match = (Match_IMP_MultiMod) matchGeneric;
			printTR();
			printTD("" + match.getPeptide().getAcidSequenceString());
			printTD("" + match.getScore());
			print ("<TD>");
			Modification [] mods = match.getModifications();
			for (int i = 0; i < mods.length; i++) {
				if (mods[i].getMonoMass() != 0) {
					print(mods[i].getMonoMass() + " at " + i + ", ");
				}
			}
			print ("</TD>");
		}
		print("</table>");
		


		printFooter();
	}
	
	
	

}
