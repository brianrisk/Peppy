package SequenceDigester;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import Peppy.Peppy;
import Peppy.NucleotideSequence;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.Sequence;
import Peppy.SequenceDigestionThread;
import Utilities.U;

public class SequenceDigester {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		U.p("Digesting sequences...");
		
		Peppy.init();
		
		//The folder where we'll store our protein files
		File digestedSequenceFolder = new File("Digested Sequences");
		digestedSequenceFolder.mkdirs();
		
		
		//Get references to our sequence files -- no nucleotide data is loaded at this point
		ArrayList<Sequence> sequences = Sequence.loadSequences(Properties.sequenceDirectoryOrFile);
		
		for (int sequenceIndex = 0; sequenceIndex < sequences.size(); sequenceIndex++)  {
			Sequence sequence = sequences.get(sequenceIndex);
			U.p("digesting sequence: " + sequence.getSequenceFile().getName());
			ArrayList<NucleotideSequence> nucleotideSequences = sequence.getNucleotideSequences();
			ArrayList<Peptide> peptides;
			for (int nucleotideSequenceIndex = 0; nucleotideSequenceIndex < nucleotideSequences.size(); nucleotideSequenceIndex++) {
				NucleotideSequence nucleotideSequence = nucleotideSequences.get(nucleotideSequenceIndex);
				for (byte frame = 0; frame < 3; frame++) {
					for (int forwards = 0; forwards <=1; forwards++) {
						SequenceDigestionThread sdt =  new SequenceDigestionThread(nucleotideSequence, frame, forwards == 0);
						peptides = sdt.getPeptides();
						File peptideFile = new File(digestedSequenceFolder, "" + sequenceIndex + "-" + nucleotideSequenceIndex + "-frame:" + frame + "-" + forwards + ".txt");
						try {
							PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(peptideFile)));
							for (Peptide peptide: peptides) {
								pw.println(peptide.getAcidSequence());
							}
							pw.flush();
							pw.close();
						} catch (IOException e) {
							e.printStackTrace();
						}
					}
				}
				
			}
		}

		U.p("done digesting.");

	}

}
