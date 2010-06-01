package Peppy;
import java.awt.Color;

/**
 * A peak is a data point in a spectrum.
 * @author Brian Risk
 *
 */
public class Peak implements Comparable<Peak>{
	
	double mass;
	double intensity;

	Color color = Color.gray;
	
	final int COMPARE_MASS = 0;
	final int COMPARE_INTENSITY = 1;
	int compareBy = COMPARE_MASS;
	
	public Peak(String s) throws MalformedPeakException {
		String [] chunks = s.split(" |\t");
		if (chunks.length != 2) {
			throw new MalformedPeakException();
		}
		mass = Double.parseDouble(chunks[0]);
		intensity = Double.parseDouble(chunks[1]);
	}
	
	public double getMass() {return mass;}
	public double getIntensity() {return intensity;}
	public Color getColor() {return color;}
	
	public void setCompareByMass() {compareBy = COMPARE_MASS;}
	public void setCompareByIntensity() {compareBy = COMPARE_INTENSITY;}
	public void setColor(Color color) {this.color = color;}
	public void setIntensity(double intensity) {this.intensity = intensity;}
	
	//@Override
	public int compareTo(Peak p) {
		if (compareBy == COMPARE_MASS) {
			if (mass < p.getMass()) return -1;
			if (mass > p.getMass()) return 1;
			return  0;
		} else {
			if (intensity < p.getIntensity()) return -1;
			if (intensity > p.getIntensity()) return 1;
			return 0;
		}
	}

}

