package Peppy;

public class Match_IMP extends Match {

	public Match_IMP(Spectrum spectrum, Peptide peptide) {
		super(spectrum, peptide);
	}

	@Override
	public void calculateScore() {
		score = -Math.log(calculateIMP());	
	}
	
	public String getScoringMethodName() {return "IMP";}
	

}
