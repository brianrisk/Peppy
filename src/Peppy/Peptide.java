package Peppy;

import Math.HasValue;


/**
 * Is a data class that stores:
 * 1) a sequence of amino acids.
 * 2) the theoretical mass
 * 3) the begin index (end index can be calculated)
 * 4) boolean forward (false = reverse)
 * 5) the reading frame
 * @author Brian Risk
 *
 */
public class Peptide implements Comparable<Peptide>, HasValue {
	
//	private String acidSequence;
	private byte [] acidSequence;
	private double mass;
	private int startIndex;
	private int stopIndex;
	private int intronStartIndex;
	private int intronStopIndex;
	private boolean forward;
	private Sequence parentSequence;
	private Protein protein;
	private boolean isSpliced;
	private boolean isMatched = false;
	private int lengthMinusOne;
	
	
	public boolean isMatched() {
		return isMatched;
	}

	public void setMatched(boolean isMatched) {
		this.isMatched = isMatched;
	}

	/**
	 * just gets an amino acid sequence.
	 * @param sequence
	 */
	public Peptide(String sequence) {
		this.acidSequence = AminoAcids.getByteArrayForString(sequence);
		this.mass = calculateMass();
		this.startIndex = 0;
		this.stopIndex = acidSequence.length;
		this.forward = true;
		this.parentSequence = null;
		this.isSpliced = false;
		this.lengthMinusOne = this.acidSequence.length - 1;
	}
	
	public Peptide(String acidSequence, int startIndex, int stopIndex, int intronStartIndex, int intronStopIndex, boolean forward, Sequence parentSequence, boolean isSpliced) {
		this.acidSequence = AminoAcids.getByteArrayForString(acidSequence);
		this.mass = calculateMass();
		this.startIndex = startIndex;
		this.stopIndex = stopIndex;
		this.intronStartIndex = intronStartIndex;
		this.intronStopIndex = intronStopIndex;
		this.forward = forward;
		this.parentSequence = parentSequence;
		this.isSpliced = isSpliced;
		this.lengthMinusOne = this.acidSequence.length - 1;
	}
	
	public Peptide(String acidSequence, int startIndex, int stopIndex, int intronStartIndex, int intronStopIndex, boolean forward, Sequence parentSequence, Protein protein, boolean isSpliced) {
		this.acidSequence = AminoAcids.getByteArrayForString(acidSequence);
		this.mass = calculateMass();
		this.startIndex = startIndex;
		this.stopIndex = stopIndex;
		this.intronStartIndex = intronStartIndex;
		this.intronStopIndex = intronStopIndex;
		this.forward = forward;
		this.parentSequence = parentSequence;
		this.protein = protein;
		this.isSpliced = isSpliced;
		this.lengthMinusOne = this.acidSequence.length - 1;
	}


	public int getLengthMinusOne() {
		return lengthMinusOne;
	}

	@Override
	public String toString() {
//		return mass + "\t" + getAcidSequenceString() + "\t" + startIndex + "\t" + proteinName;
		return  getAcidSequenceString() + "\t" + getMass() + "\t" + getStartIndex() + "\t" +  getStopIndex()  + "\t" + forward;
//		return  getAcidSequenceString();
	}


	public int compareTo(Peptide other) {
		if (mass > other.getMass()) return  1;
		if (mass < other.getMass()) return -1;
		return 0;
	}
	
	/**
	 * Okay, this equals is not in line with the way things work for compareTo.
	 * this compares acid sequences for equality.  compareTo compares masses.
	 * 
	 * the real trick for equality is ignoring any trailing stop (".") codon
	 */
	public boolean equals(byte [] otherAcidSequence) {
		int ourLength = acidSequence.length;
		int theirLength = otherAcidSequence.length;
		
		//ignoring terminating STOPs
		if (acidSequence[acidSequence.length - 1] == AminoAcids.STOP) ourLength--;
		if (otherAcidSequence[otherAcidSequence.length - 1] == AminoAcids.STOP) theirLength--;
		
		//if peptids not same length, they are not equal
		if (ourLength != theirLength) return false;
		
		//compare each acid
		for (int i = 0; i < ourLength; i++) {
			if (acidSequence[i] != otherAcidSequence[i]) return false;
		}
		
		//if reached this point, they are equal
		return true;
	}
	
	/**
	 * Equals if every sequential acid weighs the same as that of the other sequence
	 */
	public boolean equalsByAcidMasses(byte [] otherAcidSequence) {
		if (!equals(otherAcidSequence)) return false;
		
		for (int i = 0; i < acidSequence.length; i++) {
			if (AminoAcids.getWeightMono(acidSequence[i]) != AminoAcids.getWeightMono(otherAcidSequence[i])) return false;
		}
		
		return true;
			
	}
	
	public boolean equals(Peptide peptide) {
		if (mass == peptide.getMass()) {
			return equals(peptide.getAcidSequence());
		} else {
			return false;
		}
	}
	
	public boolean equals(String acidSequenceString) {
		return equals(AminoAcids.getByteArrayForString(acidSequenceString));
	}


	/**
	 * @return the sequence
	 */
	public byte [] getAcidSequence() {
		return acidSequence;
	}
	
	public String getAcidSequenceString() {
		return AminoAcids.getStringForByteArray(acidSequence);
	}
	
	public int getLength() {
		return acidSequence.length;
	}


	/**
	 * @return the mass
	 */
	public double getMass() {
		return mass;
	}


	/**
	 * @return the index
	 */
	public int getStartIndex() {
		if (forward) {
			return startIndex;
		} else {
			return stopIndex + 1;
		}
	}
	
	public int getStopIndex() {
		if (forward) {
			return stopIndex;
		} else {
			return startIndex + 1;
		}
	}

	public int getIntronStartIndex() {
		if (forward) {
			return intronStartIndex;
		} else {
			return intronStopIndex + 1;
		}
	}

	public int getIntronStopIndex() {
		if (forward) {
			return intronStopIndex;
		} else {
			return intronStartIndex + 1;
		}
	}

	public Protein getProtein() {
		return protein;
	}


//	/**
//	 * This returns the start position.  For reporting purposes.  Conforms to standards.
//	 * @return
//	 */
//	public int getSTART() {
//		if (forward) {
//			return startIndex;
//		} else {
//			return startIndex + 1 - (acidSequence.length() * 3);
//		}
//	}
//	
//	public int getSTOP() {
//		if (forward) {
//			return startIndex + (acidSequence.length() * 3);
//		} else {
//			return startIndex + 1;
//		}
//	}


	/**
	 * @return the forward
	 */
	public boolean isForward() {
		return forward;
	}


	public boolean isSpliced() {
		return isSpliced;
	}

	/**
	 * @return the parentSequence
	 */
	public Sequence getParentSequence() {
		return parentSequence;
	}
	
	/**
	 * This will calculate either the mono mass or the average mass depending on the setting
	 * in your Properties object.
	 * @return
	 */
	private double calculateMass() {
		double mass = 0.0;
		if (Properties.useMonoMass) {
			for (int i = 0; i < acidSequence.length; i++) {
				if (AminoAcids.isValid(acidSequence[i])) {
					mass += AminoAcids.getWeightMono(acidSequence[i]);
				} else {
					mass = -1;
					return mass;
				}
			}
			mass += Definitions.WATER_MONO;
		} else {
			for (int i = 0; i < acidSequence.length; i++) {
				if (AminoAcids.isValid(acidSequence[i])) {
					mass += AminoAcids.getWeightAverage(acidSequence[i]);
				} else {
					mass = -1;
					return mass;
				}
			}
			mass += Definitions.WATER_AVERAGE;
		}
		return mass;
	}

	public double getValue() {
		return getMass();
	}
	

}
