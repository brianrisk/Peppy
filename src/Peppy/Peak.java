package Peppy;
import java.awt.Color;

/**
 * A peak is a data point in a spectrum.
 * @author Brian Risk
 *
 */
public class Peak implements Comparable<Peak>{
	
	float mass;
	float intensity;
	int yIonNumber = -1;
	int bIonNumber = -1;
	
	public boolean used = false;

	Color color = Color.gray;
	boolean hilighted = false;
	
	public boolean isHilighted() {
		return hilighted;
	}

	public void setHilighted(boolean hilighted) {
		this.hilighted = hilighted;
	}

	private static int tracker = 0;
	final int COMPARE_MASS = tracker++;
	final int COMPARE_INTENSITY = tracker++;
	int compareBy = COMPARE_MASS;
	
	public Peak(String s) throws Exception {
		String [] chunks = s.split(" |\t");
		if (chunks.length != 2) {
			throw new Exception("Malformed peak exception for string: \"" + s + "\"");
		}
		mass = Float.parseFloat(chunks[0]);
		intensity = Float.parseFloat(chunks[1]);
	}
	
	public float getMass() {return mass;}
	public float getIntensity() {return intensity;}
	public Color getColor() {return color;}
	
	/**
	 * @return the yIonNumber
	 */
	public int getyIonNumber() {
		return yIonNumber;
	}

	/**
	 * @return the bIonNumber
	 */
	public int getbIonNumber() {
		return bIonNumber;
	}

	/**
	 * @param yIonNumber the yIonNumber to set
	 */
	public void setyIonNumber(int yIonNumber) {
		this.yIonNumber = yIonNumber;
	}

	/**
	 * @param bIonNumber the bIonNumber to set
	 */
	public void setbIonNumber(int bIonNumber) {
		this.bIonNumber = bIonNumber;
	}

	public void setCompareByMass() {compareBy = COMPARE_MASS;}
	public void setCompareByIntensity() {compareBy = COMPARE_INTENSITY;}
	public void setColor(Color color) {this.color = color;}
	public void setIntensity(float intensity) {this.intensity = intensity;}
	
	//@Override
	public int compareTo(Peak p) {
		if (compareBy == COMPARE_MASS) {
			//want to sort from least to greatest
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

