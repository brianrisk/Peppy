package Experimental;

import Peppy.Match;
import Peppy.MatchesSpectrum;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.Spectrum;
import Peppy.U;

public class psm {
	
	public static void main(String args[]) {
		MatchesSpectrum matchesSpectrum;
		Spectrum spectrum;
		Peptide peptide;
		Match match;
		double newMatchScore;
		String spectrumLocation = "/Users/risk2/PeppyData/ENCODE/GM12878/spectra uncompressed/wcl/15305/101111_GM_WCL_SDS_2_II.5801.5801.2.dta";
		String peptideSequence = "APAGSAAGEGLLPHR";
//		String peptideSequence = "ILPMAMLLIFSR"; 
		
		
		/* creating the PSM */
		spectrum = new Spectrum(spectrumLocation);
		matchesSpectrum = new MatchesSpectrum(spectrum);
		peptide = new Peptide(peptideSequence);
		match = Properties.matchConstructor.createMatch(matchesSpectrum, peptide);
		newMatchScore = match.getScore();
		
		for (int fragmentTolerance = 10; fragmentTolerance < 1000; fragmentTolerance += 10) {
			Properties.fragmentTolerance = fragmentTolerance;
			
			spectrum = new Spectrum(spectrumLocation);
			matchesSpectrum = new MatchesSpectrum(spectrum);
			peptide = new Peptide(peptideSequence);
			match = Properties.matchConstructor.createMatch(matchesSpectrum, peptide);

			newMatchScore = match.getScore();
			U.p(newMatchScore + 3);
		}
		
	}

}
