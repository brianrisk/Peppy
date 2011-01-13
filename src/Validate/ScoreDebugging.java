package Validate;

import java.io.File;

import Peppy.Peptide;
import Peppy.Spectrum;
import Peppy.Match;
import Utilities.U;

public class ScoreDebugging {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		U.p("debugging the score");
		
		Spectrum spectrum;
		Peptide peptide;
		Match match;
		
		
//		Spectrum spectrum = new Spectrum(new File("/Users/risk2/PeppyOverflow/tests/ecoli/spectra/C22_MSMS_1962.9415_2.txt"));
//		Peptide peptide = new Peptide("DLSELSTYSFVDNVAFR");
//		Match match = new Match(spectrum, peptide, null);
//		
//		U.p();
//		U.p("ours:");
//		spectrum = new Spectrum(new File("/Users/risk2/PeppyOverflow/tests/ecoli/spectra/C22_MSMS_1962.9415_2.txt"));
//		peptide = new Peptide("ACEDLCTIEQRPKMDGR");
//		match = new Match(spectrum, peptide, null);
		
		

	}
	
	/*
	 * also bad.  it looks like every b ion is also a y.  A symmetrical fucking theoretical spectrum.
	 * E value: 1.8080197160165963 
Spectrum file name: E10_MSMS_1066.6492_6.t2d.txt 
Correct peptide is in the database: true 
our score: 12, 25.426829064993225 
their score: 11, 24.94633213500149 
ours vs theirs: SSKSSNSSKR vs MKPFIFGAR
	 */

}
