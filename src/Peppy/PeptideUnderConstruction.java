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
	
	public PeptideUnderConstruction(int startIndex, char aminoAcid) {
		this.startIndex = startIndex;
		addAminoAcid(aminoAcid);
	}
	
	public PeptideUnderConstruction(int sequenceIndex) {
		this.startIndex = sequenceIndex;
	}
	
	public void addAminoAcid(char acid) {
		buffer.append(acid);
		if (acid == '.') breakCount++;
		if (acid == 'K') breakCount++;
		if (acid == 'R') breakCount++;
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

}
