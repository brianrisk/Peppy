package Peppy;

/**
 * Used for loading in Mascot and X!Tandem results
 * 
 * We already have a score, E value, etc for a match
 * and simply need to contain it in an object
 * 
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class Match_Blank extends Match{
	
	public Match_Blank(Spectrum spectrum, Peptide peptide, double score) {
		this.matchesSpectrum = new MatchesSpectrum(spectrum);
		this.peptide = peptide;
		this.score = score;
	}

	@Override
	public void calculateScore() {
		//do nothing
	}


}
