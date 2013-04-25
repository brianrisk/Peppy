package Peppy;


/**
 * Copyright 2012, Brian Risk
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
