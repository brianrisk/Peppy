package Peppy;

/**
 * A peak is a data point in a spectrum.
 * 
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class Peak implements Comparable<Peak>{
	
	float mass;
	float intensity;
	
	public boolean used = false;
	
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
	
	public Peak(float mass, float intensity) {
		this.mass = mass;
		this.intensity = intensity;
	}
	
	public float getMass() {return mass;}
	
	public float getIntensity() {return intensity;}

	public void setCompareByMass() {compareBy = COMPARE_MASS;}
	
	public void setCompareByIntensity() {compareBy = COMPARE_INTENSITY;}
	
	public void setIntensity(float intensity) {this.intensity = intensity;}
	
	public void setMass(float mass) {this.mass = mass;}
	
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

