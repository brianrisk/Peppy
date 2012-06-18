package Math;

/**
 * Mass error can be measured absolutely in Daltons (atomic mass units)
 * or in PPM (parts per million).  PPM is based on the mass being measured.
 * @author Brian Risk
 *
 */
public class MassError {
	
	
	/**
	 * based on the mass and PPM, this finds our error in Daltons
	 * @param mass
	 */
	public static double getDaltonError(double PPM, double mass) {
		return mass * PPM / 1000000;
//		return 1000 * PPM / 1000000; 
	}
	
	
	/**
	 * If we want to see how many PPM difference from a theoretical peptide (true) and
	 * a spectrum (observed)
	 * @return
	 */
	public static double getPPMDifference(double trueValue, double observedValue) {
		return ((trueValue - observedValue) / trueValue) * 1000000;
	}
	


}
