package Peppy;

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
public class Peptide implements Comparable<Peptide> {
	
	private String acidSequence;
	private double mass;
	private int startIndex;
	private int stopIndex;
	private int intronStartIndex;
	private int intronStopIndex;
	private boolean forward;
	private byte readingFrame;
	private Sequence parentSequence;
	private String proteinName;
	private boolean isSpliced;
	
	
	/**
	 * just gets an amino acid sequence.
	 * @param sequence
	 */
	public Peptide(String sequence) {
		this.acidSequence = sequence;
		this.mass = calculateMass();
		this.startIndex = 0;
		this.forward = true;
		this.readingFrame = (byte) 0;
		this.parentSequence = null;
		this.isSpliced = false;
	}
	
	/**
	 * A version of the constructor which calculates the mass from the given sequence.
	 * @param acidSequence
	 * @param startIndex
	 * @param forward
	 */
	public Peptide(String acidSequence, int startIndex, boolean forward, byte readingFrame, Sequence parentSequence) {
		this.acidSequence = acidSequence;
		this.mass = calculateMass();
		this.startIndex = startIndex;
		this.forward = forward;
		this.readingFrame = readingFrame;
		this.parentSequence = parentSequence;
		this.isSpliced = false;
	}
	
	/**
	 * For spliced peptides
	 * @param acidSequence
	 * @param startIndex
	 * @param forward
	 */
	public Peptide(String acidSequence, int startIndex, int stopIndex, boolean forward, byte readingFrame, Sequence parentSequence, boolean isSpliced) {
		this.acidSequence = acidSequence;
		this.mass = calculateMass();
		this.startIndex = startIndex;
		this.stopIndex = stopIndex;
		this.forward = forward;
		this.readingFrame = readingFrame;
		this.parentSequence = parentSequence;
		this.isSpliced = isSpliced;
	}
	
	/**
	 * A constructor for if our peptide comes from a protein database, not from digested DNA
	 * @param acidSequence
	 * @param proteinName
	 */
	public Peptide(String acidSequence, String proteinName) {
		this.acidSequence = acidSequence;
		this.mass = calculateMass();
		this.proteinName = proteinName;
	}


	public String getProteinName() {
		return proteinName;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
//		return "forward=" + forward + ", index=" + index + ", mass="
//				+ mass + ", readingFrame=" + readingFrame + ", sequence="
//				+ sequence;
		int outFrame = readingFrame + 1;
		if (!forward) outFrame *= -1;
		return mass + "\t" + acidSequence + "\t" + startIndex + "\t" + outFrame + "\t" + proteinName;
	}


	public int compareTo(Peptide o) {
		if (mass > o.getMass()) return 1;
		if (mass < o.getMass()) return -1;
		return 0;
	}
	
	/**
	 * Okay, this equals is not in line with the way things work for compareTo.
	 * this compares acid sequences for equality.  compareTo compares masses.
	 * 
	 * the real trick for equality is ignoring any trailing stop (".") codon
	 */
	public boolean equals(String otherAcidSequence) {
		String us = acidSequence.toUpperCase();
		if (acidSequence.endsWith(".")) {
			us = us.substring(0, us.indexOf('.'));
		}
		
		String them = otherAcidSequence.toUpperCase();
		if (otherAcidSequence.endsWith(".")) {
			them = them.substring(0, them.indexOf('.'));
		}
		
		if (us.length() != them.length()) return false;
		
		return us.equals(them);
			
	}
	
	/**
	 * Equals if, by mass, each amino acid is equal
	 * 
	 * the real trick for equality is ignoring any trailing stop (".") codon
	 */
	public boolean equalsByAcidMasses(String otherAcidSequence) {
		String us = acidSequence.toUpperCase();
		if (acidSequence.endsWith(".")) {
			us = us.substring(0, us.indexOf('.'));
		}
		
		String them = otherAcidSequence.toUpperCase();
		if (otherAcidSequence.endsWith(".")) {
			them = them.substring(0, them.indexOf('.'));
		}
		
		if (us.length() != them.length()) return false;
		
		//go through each amino acid and, if they weigh the same, they are considered equal
		boolean equal = true;
		
		for (int i = 0; i < us.length(); i++) {
			if (Definitions.getAminoAcidWeightMono(us.charAt(i)) != Definitions.getAminoAcidWeightMono(them.charAt(i))) {
				equal = false;
				break;
			}
		}
		return equal;
			
	}
	
	public boolean equals(Peptide peptide) {
		if (mass == peptide.getMass()) {
			return equals(peptide.getAcidSequence());
		} else {
			return false;
		}
	}


	/**
	 * @return the sequence
	 */
	public String getAcidSequence() {
		return acidSequence;
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
		return startIndex;
	}
	
	public int getStopIndex() {
		return stopIndex;
	}

	public int getIntronStartIndex() {
		return intronStartIndex;
	}

	public int getIntronStopIndex() {
		return intronStopIndex;
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
	 * @return the readingFrame
	 */
	public byte getReadingFrame() {
		return readingFrame;
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
	public double calculateMass() {
		double mass = 0.0;
		if (Properties.useMonoMass) {
			for (int i = 0; i < acidSequence.length(); i++) {
				mass += Definitions.getAminoAcidWeightMono(acidSequence.charAt(i));
			}
			mass += Definitions.WATER_MONO;
		} else {
			for (int i = 0; i < acidSequence.length(); i++) {
				mass += Definitions.getAminoAcidWeightAverage(acidSequence.charAt(i));
			}
			mass += Definitions.WATER_AVERAGE;
		}
		return mass;
	}

}
