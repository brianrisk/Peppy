package Peppy;

import Utilities.U;

public class MatchConstructor {
	
	Class<? extends Match> match;
	
	@SuppressWarnings("unchecked")
	public MatchConstructor(String className) {
		try {
			match = (Class<? extends Match>) Class.forName(className);
		} catch (ClassNotFoundException e) {
			U.p("A match classname was not recognized, yo");
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	public Match createMatch(Spectrum spectrum, Peptide peptide) {
		return createMatch(spectrum, peptide, true);
	}
	
	public Match createMatch(Spectrum spectrum, Peptide peptide, boolean calculateScore) {
		Match out = null;
		try {
			out = match.newInstance();
			out.setSpectrum(spectrum);
			out.setPeptide(peptide);
			if (calculateScore) out.calculateScore();
		} catch (InstantiationException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return out;
	}

}
