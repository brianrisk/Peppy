package Peppy;

import java.util.ArrayList;

/**
 * A data structure to hold a string of DNA nucleotides
 * 
 * Copyright 2012, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class NucleotideSequence {
	
	private String sequenceDescription;
	private String sequence;
	private Sequence_DNA parentSequence;
	
	public static void main (String args[]) {
		NucleotideSequence ns = new NucleotideSequence(">test", "ATTTGTCTTCGATGACATCAACAAGAGCAAGTTCATCTGCCAAGGC", null);
		ArrayList<NucleotideSequence> array = new ArrayList<NucleotideSequence>();
		array.add(ns);
//		for (int i = 0; i < 12; i++) {array.add(ns);}
		ShortNucleotideDigestionServer server = new ShortNucleotideDigestionServer(array, false);
		ArrayList<Protein> proteins = server.getResults();
		U.p(proteins.size());
		for (Protein protein: proteins) {
			U.p(protein.getAcidString());
		}
	}
	
	public NucleotideSequence(String sequenceDescription, String sequence, Sequence_DNA parentSequence) {
		this.sequenceDescription = sequenceDescription;
		this.sequence = sequence;
		this.parentSequence = parentSequence;
	}
	
	public String getChromosomeNumberString() {
		return sequenceDescription.substring(">chr".length());
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
