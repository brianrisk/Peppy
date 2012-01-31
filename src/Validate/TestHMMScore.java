package Validate;

import Peppy.Match;
import Peppy.MatchConstructor;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.Spectrum;
import Peppy.U;

public class TestHMMScore {
	
	public static void main(String args[]) {
		Properties.scoringMethodName = "Peppy.Match_HMM";
		Properties.calculateEValues = true;
		Properties.highIntensityCleaning = true;
		Properties.matchConstructor = new MatchConstructor(Properties.scoringMethodName);
		
		Spectrum spectrum;
		Peptide peptide;
		Match match;
		
		/* a great match */
		spectrum = new Spectrum("/Users/risk2/PeppyOverflow/tests/human/spectra/JNA-and-JGP-AAP-27Aug03_HUPO.1776.1776.2.dta");
		peptide = new Peptide("LYHSEAFTVNFGDTEEAKK");
		match = Properties.matchConstructor.createMatch(spectrum, peptide);
		U.p("IMP score is: " + match.getIMP());
		U.p("HMM score is: " + match.getScore());
		
		
		/* a not great match */
		spectrum = new Spectrum("/Users/risk2/PeppyOverflow/tests/aurum/spectra/T10707_Well_F23_1841.99_19185.mgf..pkl");
		peptide = new Peptide("VLPFLLEFWQGQPNR");
		match = Properties.matchConstructor.createMatch(spectrum, peptide);
		U.p("IMP score is: " + match.getIMP());
		U.p("HMM score is: " + match.getScore());
		
		/* a reversed match */
		spectrum = new Spectrum("/Users/risk2/PeppyOverflow/tests/aurum/spectra/T10707_Well_F23_1841.99_19185.mgf..pkl");
		peptide = new Peptide("RNPQGQWFELLFPLV");
		match = Properties.matchConstructor.createMatch(spectrum, peptide);
		U.p("IMP score is: " + match.getIMP());
		U.p("HMM score is: " + match.getScore());
		
		/* a reversed match */
		spectrum = new Spectrum("/Users/risk2/PeppyOverflow/tests/aurum/spectra/T10707_Well_F23_1841.99_19185.mgf..pkl");
		peptide = new Peptide("RFELPQLWGNLQFPV");
		match = Properties.matchConstructor.createMatch(spectrum, peptide);
		U.p("IMP score is: " + match.getIMP());
		U.p("HMM score is: " + match.getScore());
		
	}

}
