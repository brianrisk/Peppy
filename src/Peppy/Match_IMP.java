package Peppy;


/**
 * Copyright 2012, Brian Risk
 * Released under the Netscape Public License
 * 
 * @author Brian Risk
 *
 */
public class Match_IMP extends Match {


	@Override
	public void calculateScore() {
		score = -Math.log10(calculateIMP());	
	}
	
	public String getScoringMethodName() {return "IMP";}


}
