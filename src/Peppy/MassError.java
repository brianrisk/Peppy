package Peppy;

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
	}
	


}
