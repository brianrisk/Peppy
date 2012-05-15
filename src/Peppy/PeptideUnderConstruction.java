package Peppy;

/**
 * This class is basically a tool to help simplify the whole
 * nucleotide translation / digestion process.
 * @author Brian Risk
 *
 */
public class PeptideUnderConstruction {
	
	private int breakCount = 0;
	private int startIndex;
	private StringBuffer buffer = new StringBuffer();
	private boolean inORF;
	private int ORFSize;
	
	public PeptideUnderConstruction(int startIndex, char aminoAcid, boolean inORF, int ORFSize) {
		this.startIndex = startIndex;
		addAminoAcid(aminoAcid);
		this.inORF = inORF;
		this.ORFSize = ORFSize;
	}
	
	public PeptideUnderConstruction(int sequenceIndex) {
		this.startIndex = sequenceIndex;
	}
	
	public void addAminoAcid(char acid) {
		buffer.append(acid);
		if (Protein.isBreak(acid)) breakCount++;
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

}
