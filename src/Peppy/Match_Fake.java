package Peppy;

import java.util.Random;

public class Match_Fake extends Match {
	
	private static Random random = new Random();

	@Override
	public void calculateScore() {
		byte [] acidSequence = peptide.getAcidSequence();
		if (acidSequence[acidSequence.length - 1] == AminoAcids.STOP && random.nextInt(30) == 1) {
			score = Math.abs(random.nextGaussian()) + 50;
		} else {
			score = Math.abs((random.nextGaussian()) * 4 + 5);
//			if (score > 15) score = 0;

		}
		
	}

	@Override
	public String getScoringMethodName() {
		return "Match_Fake";
	}

}
