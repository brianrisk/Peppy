package Peppy;

/**
 * This class is for when we already have a score, E value, etc for a match
 * and simply need to contain it in an object
 * @author Brian Risk
 *
 */
public class Match_Blank extends Match{
	
	public Match_Blank(Spectrum spectrum, Peptide peptide, double score, double eValue) {
		this.spectrum = spectrum;
		this.peptide = peptide;
		this.score = score;
		this.eValue = eValue;
		this.rank = 1;
	}

	@Override
	public void calculateScore() {
		//do nothing
	}

	@Override
	public String getScoringMethodName() {
		return "blank";
	}

}
