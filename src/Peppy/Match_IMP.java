package Peppy;

public class Match_IMP extends Match {


	@Override
	public void calculateScore() {
		score = -Math.log(calculateIMP());	
	}
	
	public String getScoringMethodName() {return "IMP";}
	

}