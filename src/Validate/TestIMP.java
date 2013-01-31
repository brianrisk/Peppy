package Validate;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import Peppy.Match;
import Peppy.MatchConstructor;
import Peppy.MatchesSpectrum;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.Spectrum;
import Peppy.U;

public class TestIMP {
	
	public static void main(String args[]) {
		Peppy.Peppy.init(args);
		
		Properties.precursorTolerance = 2000;
		Properties.fragmentTolerance = 300;
		Properties.iodoacetamideDerivative = true;
		
		/* What scoring mechanism? */
		Properties.scoringMethodName = "Peppy.Match_IMP";
		Properties.matchConstructor = new MatchConstructor(Properties.scoringMethodName);
		
		testIndividuals();

	}
	
	/**
	 * This takes a Validation report test file and sees if the scoring system has improved at all
	 */
	public static void testFile() {
		
		MatchesSpectrum matchesSpectrum;
		Spectrum spectrum;
		Peptide peptide;
		Match match;
		
		String spectrumBase = "/Users/risk2/PeppyData/tests/USP top 10/spectra/";
		
		File file = new File("/Users/risk2/Documents/workspace/JavaGFS/validation/24 1359399551658/USP top 10/USP top 10 results.txt");
		
		int oldTrueCount = 0;
		int newTrueCount = 0;
		int lineCount = 0;
		final double fdr = 0.01;
		double oldOnePercent = 0;
		double newOnePercent = 0;
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			String line = br.readLine();
			line = br.readLine();
			
			
			while (line != null) {
				String [] chunks = line.split("\t");
				/*
				 * index 0: true match or not
				 * index 1: id
				 * index 2: match score
				 * index 3: match peptide
				 * index 4: true score
				 * index 5: true peptide
				 * index 6: spectrum name
				 */
				
				if (chunks[0].equals("true")) oldTrueCount++;
				
				spectrum = new Spectrum(spectrumBase + chunks[6]);
				matchesSpectrum = new MatchesSpectrum(spectrum);
				peptide = new Peptide(chunks[3]);
				match = Properties.matchConstructor.createMatch(matchesSpectrum, peptide);
				double newMatchScore = match.getScore();
				
				spectrum = new Spectrum(spectrumBase + chunks[6]);
				matchesSpectrum = new MatchesSpectrum(spectrum);
				peptide = new Peptide(chunks[5]);
				match = Properties.matchConstructor.createMatch(matchesSpectrum, peptide);
				double newTrueScore = match.getScore();
				
				if (newTrueScore >= newMatchScore) newTrueCount++;
				
				/* calculate FDR */ 
				lineCount++;
				if ((double) (lineCount - oldTrueCount) / lineCount < fdr) {
					oldOnePercent = lineCount;
				}
				if ((double) (lineCount - newTrueCount) / lineCount < fdr) {
					newOnePercent = lineCount;
				}
				
				line = br.readLine();
				
			}
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		if (newTrueCount > oldTrueCount) {
			U.p("improvement!");
		}
		U.p("oldTrueCount: " + oldTrueCount);
		U.p("newTrueCount: " + newTrueCount);
		U.p("oldFDR: " + (oldOnePercent  * 100 / lineCount));
		U.p("newFDR: " + (newOnePercent  * 100 / lineCount));
		
	}
	
	public static void testIndividuals() {
		MatchesSpectrum matchesSpectrum;
		Spectrum spectrum;
		Peptide peptide;
		Match match;
		
		/* a great match */
		spectrum = new Spectrum("/Users/risk2/PeppyData/tests/USP top 10/spectra/USP_std_run1_100504135936.7340.7340.2.dta");
		matchesSpectrum = new MatchesSpectrum(spectrum);
		peptide = new Peptide("AFLEVNEEGSEAAASTAVVIAGR");
		match = Properties.matchConstructor.createMatch(matchesSpectrum, peptide);
		U.p("score is: " + match.getScore());
		
		
//		spectrum = new Spectrum("/Users/risk2/PeppyData/tests/USP top 10/spectra/USP_std_run1_100504135936.5338.5338.2.dta");
//		matchesSpectrum = new MatchesSpectrum(spectrum);
//		peptide = new Peptide("YYTLEEIKK");
//		match = Properties.matchConstructor.createMatch(matchesSpectrum, peptide);
//		U.p("score is: " + match.getScore());
//		
//		spectrum = new Spectrum("/Users/risk2/PeppyData/tests/USP top 10/spectra/USP_std_run1_100504135936.5338.5338.2.dta");
//		matchesSpectrum = new MatchesSpectrum(spectrum);
//		peptide = new Peptide("YYTLEEIQK");
//		match = Properties.matchConstructor.createMatch(matchesSpectrum, peptide);
//		U.p("score is: " + match.getScore());
//		
//		U.p();
//		
//		spectrum = new Spectrum("/Users/risk2/PeppyData/tests/USP top 10/spectra/USP_std_run1_100504135936.4792.4792.2.dta");
//		matchesSpectrum = new MatchesSpectrum(spectrum);
//		peptide = new Peptide("FEELLTR");
//		match = Properties.matchConstructor.createMatch(matchesSpectrum, peptide);
//		U.p("score is: " + match.getScore());
//		
//		spectrum = new Spectrum("/Users/risk2/PeppyData/tests/USP top 10/spectra/USP_std_run1_100504135936.4792.4792.2.dta");
//		matchesSpectrum = new MatchesSpectrum(spectrum);
//		peptide = new Peptide("SLGVGFATR");
//		match = Properties.matchConstructor.createMatch(matchesSpectrum, peptide);
//		U.p("score is: " + match.getScore());

	}

}
