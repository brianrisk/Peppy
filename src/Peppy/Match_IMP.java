package Peppy;


/**
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class Match_IMP extends Match {


	@Override
	public void calculateScore() {
		score = -Math.log10(calculateIMP());	
	}
	


}
