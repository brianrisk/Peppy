package Peppy;

/**
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class PeakLabeled extends Peak {
	
	String color = "#000000";
	String label = "";
	int bIonNumber = 0;
	int yIonNumber = 0;

	public PeakLabeled(Peak peak) {
		super(peak.getMass(), peak.getIntensity());
		// TODO Auto-generated constructor stub
	}

	public String getColor() {
		return color;
	}

	public void setColor(String color) {
		this.color = color;
	}

	public String getLabel() {
		return label;
	}

	public void setLabel(String label) {
		this.label = label;
	}

	public int getbIonNumber() {
		return bIonNumber;
	}

	public void setbIonNumber(int bIonNumber) {
		this.bIonNumber = bIonNumber;
	}

	public int getyIonNumber() {
		return yIonNumber;
	}

	public void setyIonNumber(int yIonNumber) {
		this.yIonNumber = yIonNumber;
	}
	
	
	public void setMass(float mass) {
		this.mass = mass;
	}

}
