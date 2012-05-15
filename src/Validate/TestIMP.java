package Validate;

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
		Properties.scoringMethodName = "Peppy.Match_IMP_VariMod";
		Properties.matchConstructor = new MatchConstructor(Properties.scoringMethodName);
		
		MatchesSpectrum matchesSpectrum;
		Spectrum spectrum;
		Peptide peptide;
		Match match;
		
		/* a great match */
		spectrum = new Spectrum("//Users/risk2/PeppyData/WashU/spectra/WHIM16/Ellis_033_WHIM2_WHIM16-frac_dta/Ellis_033_2700_261_11/Ellis_033_2700_261_11.mgf.12910.dta");
		matchesSpectrum = new MatchesSpectrum(spectrum);
		peptide = new Peptide("MDEAYMNK");
		match = Properties.matchConstructor.createMatch(matchesSpectrum, peptide);
		U.p("has modification: " + match.hasModification());
		U.p("modification Mass: " + match.getMoificationdMass());
		U.p("score is: " + match.getScore());
		
		U.p();
		
		spectrum = new Spectrum("//Users/risk2/PeppyData/WashU/spectra/WHIM16/Ellis_033_WHIM2_WHIM16-frac_dta/Ellis_033_2700_261_11/Ellis_033_2700_261_11.mgf.12910.dta");
		matchesSpectrum = new MatchesSpectrum(spectrum);
		peptide = new Peptide("DVDEAYMNK");
		match = Properties.matchConstructor.createMatch(matchesSpectrum, peptide);
		U.p("has modification: " + match.hasModification());
		U.p("modification Mass: " + match.getMoificationdMass());
		U.p("score is: " + match.getScore());
		
		
		
		

		
		
	}

}
