package Peppy;


/**
 * 
 * I'm setting up "squared" to be a way to help sort out the situation where 
 * matching 20 out of 50 AAs in a peptide produces a better score than matching
 * 10 out of 10 AAs in a shorter peptide.
 * 
 * One may think percent matched would do the trick, but that would overly
 * favor short peptides.  Then a PSM that matched the 8 ions of a 4 AA peptide
 * would outscore one that matched 59 ions of a 30 AA peptide.
 * 
 * This could, potentially, be a method to eliminate the false positives produced
 * by modified matches outscoring the modified that have equal quality of match
 * 
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */

public class Match_Squared extends Match {
	
	

	@Override
	public void calculateScore() {
		int numberOfIonsMatched = getNumberOfIonsMatched();
		score = numberOfIonsMatched * numberOfIonsMatched / peptide.getLengthMinusOne();
	}


}
