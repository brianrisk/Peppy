package Peppy;

/**
 * A data structure to hold a string of DNA nucleotides
 * @author Brian Risk
 *
 */
public class DNA_Sequence {
	
	private String sequenceDescription;
	private String sequence;
	private Sequence parentSequence;
	

	
	public DNA_Sequence(String sequenceDescription, String sequence, Sequence parentSequence) {
		this.sequenceDescription = sequenceDescription;
		this.sequence = sequence;
		this.parentSequence = parentSequence;
	}
	

	public String getSequenceDescription() {
		return sequenceDescription;
	}

	public String getSequence() {
		return sequence;
	}
	
	public Sequence getParentSequence() {
		return parentSequence;
	}

}
