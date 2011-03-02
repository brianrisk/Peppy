package Peppy;

import HMMScore.HMMClass;

public class Match_HMM extends Match {

	public Match_HMM(Spectrum spectrum, Peptide peptide) {
		super(spectrum, peptide);
	}

	@Override
	public void calculateScore() {
		HMMClass scorer = new HMMClass(peptide.getAcidSequenceString(), spectrum);
		score = scorer.score();
	}
	
	public String getScoringMethodName() {return "HMM_score";}


}
