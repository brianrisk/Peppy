package USP;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

import Peppy.JavaGFS;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.ProteinDigestion;
import Peppy.Sequence;
import Peppy.Spectrum;
import Peppy.Match;
import Utilities.U;
import Validate.ReliabilityTester;

public class TuneMSMSFit {
	
	public static void main(String args[]) {
		//print intro
		U.p("Hello!  We're tuning up MSMS Fit for the USP data!");
		
		//load the correct peptide set
		ArrayList<Peptide> correctPeptides = ProteinDigestion.getPeptidesFromProteinFile(new File("USP/extracted-proteins.txt"));
		
		//load the high scoring peptide set
		ArrayList<Peptide> peptides = ReliabilityTester.loadHighScoringPeptides("USP");
		
		//load our spectra
		ArrayList<Spectrum> spectra = Spectrum.loadSpectraFromFolder("tests/USP/spectra");
		
		//Go through each spectra and find the best matches
		Properties.peakDifferenceThreshold = 0.25;
		Properties.maximumNumberOfMatchesForASpectrum = 1;
		ArrayList<Match> matches = JavaGFS.asynchronousDigestion(peptides, spectra, null);
		
		//go through each of the best matches and find how many of them exist in the correct peptide set
		int hits = 0;
		for (Match match: matches) {
			for (Peptide peptide: correctPeptides) {
				if (match.getPeptide().getAcidSequence().equals(peptide.getAcidSequence())) {
					hits++;
					break;
				}
			}
		}
		
		//print results
		U.p("we got this many hits: " + hits);
		U.p("that's this percent: " + (double) hits / spectra.size() * 100);
	}
	
	
}
