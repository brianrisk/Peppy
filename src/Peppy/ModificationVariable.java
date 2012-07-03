package Peppy;

/**
 * For the purpose 
 * 
 * Copyright 2012, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class ModificationVariable {
	
	byte AminoAcid;
	double mass;
	
	public ModificationVariable(byte aminoAcid, double mass) {
		AminoAcid = aminoAcid;
		this.mass = mass;
	}

	public byte getAminoAcid() {
		return AminoAcid;
	}

	public double getMass() {
		return mass;
	}
	

}
