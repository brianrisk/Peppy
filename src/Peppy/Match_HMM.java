package Peppy;

import HMMScore.HMMClass;

public class Match_HMM extends Match {
	
	static {
		HMMScore.HMMClass.HmmSetUp();
	}

	@Override
	public void calculateScore() {
		HMMClass scorer = new HMMClass(peptide.getAcidSequenceString(), spectrum);
		score = scorer.score();
	}
	
	@Override
	public String getScoringMethodName() {return "HMM_score";}


}
