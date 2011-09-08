package Peppy;

/**
 * A data structure to hold a string of DNA nucleotides
 * @author Brian Risk
 *
 */
public class Nucleotides {
	
	private String sequenceDescription;
	private String sequence;
	private Sequence_DNA parentSequence;
	

	
	public Nucleotides(String sequenceDescription, String sequence, Sequence_DNA parentSequence) {
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
	
	public Sequence_DNA getParentSequence() {
		return parentSequence;
	}
	
	public void clearSequenceData() {
		sequence = null;
	}

}
