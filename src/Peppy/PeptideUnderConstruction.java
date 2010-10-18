package Peppy;

/**
 * This class is basically a tool to help simplify the whole
 * nucleotide translation / digestion process.
 * @author Brian Risk
 *
 */
public class PeptideUnderConstruction {
	
	private int breakCount = 0;
	private int sequenceIndex;
	private StringBuffer buffer = new StringBuffer();
	
	public PeptideUnderConstruction(int sequenceIndex, char aminoAcid) {
		this.sequenceIndex = sequenceIndex;
		addAminoAcid(aminoAcid);
	}
	
	public PeptideUnderConstruction(int sequenceIndex) {
		this.sequenceIndex = sequenceIndex;
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
	
	public int getCodeChunkIndex() {
		return sequenceIndex;
	}

}
