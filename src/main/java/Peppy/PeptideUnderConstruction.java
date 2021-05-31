package Peppy;

/**
 * This class is basically a tool to help simplify the whole
 * nucleotide translation / digestion process.
 * 
 * Copyright 2013, Brian Risk
 * 
 * 
 * 
 * @author Brian Risk
 *
 */
public class PeptideUnderConstruction {
	
	private int breakCount = 0;
	private int startIndex;
	private StringBuffer buffer = new StringBuffer();
	private boolean inORF;
	private int ORFSize;
	private char previousAminoAcid;
	
	public PeptideUnderConstruction(int startIndex, char aminoAcid, char nextAminoAcid, boolean inORF, int ORFSize, char previousAminoAcid) {
		this.startIndex = startIndex;
		addAminoAcid(aminoAcid, nextAminoAcid);
		this.inORF = inORF;
		this.ORFSize = ORFSize;
		this.previousAminoAcid = previousAminoAcid;
	}
	
	public PeptideUnderConstruction(int sequenceIndex) {
		this.startIndex = sequenceIndex;
	}
	
	public void addAminoAcid(char acid, char nextAminoAcid) {
		buffer.append(acid);
		if (Protein.isBreak(acid, nextAminoAcid)) breakCount++;
	}
	
	public int getBreakCount() {
		return breakCount;
	}
	
	public String getSequence() {
		return buffer.toString();
	}
	
	public int getStartIndex() {
		return startIndex;
	}

	public boolean isInORF() {
		return inORF;
	}
	
	public int getORFSize() {
		return ORFSize;
	}
	
	public char getPreviousAminoAcid() {
		return previousAminoAcid;
	}

}
