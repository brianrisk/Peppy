package Peppy;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

import Utilities.U;

/**
 * A "sequence" is a DNA sequence, typically in a FASTA file
 * This farms out the digestion of the nucleotide sequences to
 * SequenceDigestionTread objects
 * @author Brian Risk
 *
 */
public class Sequence_DNA extends Sequence{

	private ArrayList<DNA_Sequence> nucleotideSequences = null;
	private int startIndex = 0;
	private int stopIndex = 0;
	private int nucleotideSequenceIndex = 0;
	//To avoid missing some peptides when digestion in discrete windows we
	//are allowing for some overlap
	private final int digestionFrameOverlap = 999; //divisible by three
	
	public Sequence_DNA(File sequenceFile) {
		this.sequenceFile = sequenceFile;
	}
	
	public Sequence_DNA(String sequenceFileName) {
		this(new File(sequenceFileName));
	}
	
	/**
	 * The amount of peptides produced by a sequence may be huge.
	 * therefore, this method is set up so that it should be called
	 * multiple times, each time it is called it returns another portion of the
	 * database of peptides.  When there are no more peptides from that sequence
	 * then null is returned.
	 * @param reverse This is if we want to reverse our database to get a null database
	 * @return a sorted ArrayList of amino acid sequence fragments from the given sequence file
	 */
	public ArrayList<Peptide> extractMorePeptides(boolean reverse) {
		//if we have yet to read in all the ATGC data, do that now
		if (nucleotideSequences == null) getNucleotideSequences();
		
		//Get whatever sequence we're working on.  There will most often only be
		//one sequence in nucleotideSequences, but sometimes more
		DNA_Sequence nucleotideSequence = nucleotideSequences.get(nucleotideSequenceIndex);
		
		//if this is the first time we've tried to extract peptides
		if (!Properties.useSequenceRegion) {
			if (stopIndex == 0) {
				stopIndex = startIndex + Properties.digestionWindowSize;
			} else {
			
				//if we've reached the end of the nucleotide sequence
				if (stopIndex == nucleotideSequence.getSequence().length()) {
					startIndex = 0;
					stopIndex = Properties.digestionWindowSize;
					nucleotideSequenceIndex++;
					if (nucleotideSequenceIndex >= nucleotideSequences.size()) {
						//returning null as signal that we have reached the end
						return null;
					} else {
						nucleotideSequence = nucleotideSequences.get(nucleotideSequenceIndex);
					}
				} else {
					startIndex += Properties.digestionWindowSize;
					stopIndex += Properties.digestionWindowSize;
				}
			}
			//overshoot correction
			if (stopIndex > nucleotideSequence.getSequence().length()) {
				stopIndex = nucleotideSequence.getSequence().length();
			}
		} else {
			//if we've already extracted the region
			if (stopIndex == Properties.sequenceRegionStop) {
				//signal that we've finished
				return null;
			} else {
				startIndex = Properties.sequenceRegionStart;
				stopIndex = Properties.sequenceRegionStop;
			}
		}
		
		
		ArrayList<Peptide> peptides = new ArrayList<Peptide>();
		
		//Create our SequenceDigestionThread ArrayList
		ArrayList<DNA_DigestionThread> digestors = new ArrayList<DNA_DigestionThread>();
		for (byte frame = 0; frame < 3; frame++) {
			digestors.add(new DNA_DigestionThread(nucleotideSequence, frame, true, startIndex - digestionFrameOverlap, stopIndex, reverse));
			if (!Properties.useOnlyForwardsFrames) {
				digestors.add(new DNA_DigestionThread(nucleotideSequence, frame, false, startIndex - digestionFrameOverlap, stopIndex, reverse));
			}
		}
		
		//create the threads and start them engines!
		ArrayList<Thread> threads = new ArrayList<Thread>(); 
		for (DNA_DigestionThread digestor: digestors) {
			Thread thread = new Thread(digestor);
			thread.start();
			threads.add(thread);
		}
		
		//Wait for them all to finish
		for (Thread thread: threads) {
			try {
				thread.join();
			} catch (InterruptedException e) {
				U.p("Digestion thread interrupted!  Bad!");
				e.printStackTrace();
			}
		}
		
		//harvest all digested peptides
		for (int digestorIndex = 0; digestorIndex < digestors.size(); digestorIndex++) {
			DNA_DigestionThread digestor = digestors.get(digestorIndex);
			peptides.addAll(digestor.getPeptides());
		}
			
		Collections.sort(peptides);
		return peptides;
	}
	
	/**
	 * Runs extractMorePeptides() until there are no more peptides
	 * to extract.  Returns the full list.
	 * @param reverse This is if we want to reverse our database to get a null database
	 * @return
	 */
	public ArrayList<Peptide> extractAllPeptides(boolean reverse) {
		ArrayList<Peptide> peptides = new ArrayList<Peptide>();
		ArrayList<Peptide> peptidesToAdd = extractMorePeptides(reverse);
		while (peptidesToAdd != null) {
			peptides.addAll(peptidesToAdd);
			peptidesToAdd = extractMorePeptides(reverse);
		}
		return peptides;
	}
	
	
	/**
	 * in one FASTA file there may be many sequences.
	 * this method returns a ArrayList containing all of the sequences
	 * @return
	 */
	public ArrayList<DNA_Sequence> getNucleotideSequences() {
		if (nucleotideSequences == null) {
			nucleotideSequences = new ArrayList<DNA_Sequence>();
			try {
				BufferedReader br = new BufferedReader(new FileReader(sequenceFile));
				String line = br.readLine();
				//lines beginning with ">" name the sequence
				String sequenceDescription = null;
				StringBuffer sequence = new StringBuffer();
				while (line != null) {
					//lines that begin with ">" are comments which name the sequence
					if (line.startsWith(">")) {
						//if sequenceDescription has not been defined, that means it is the first in the file
						if (sequenceDescription != null) {
							nucleotideSequences.add(new DNA_Sequence(sequenceDescription, sequence.toString().toUpperCase(), this));
							sequence = new StringBuffer();
						}
						sequenceDescription = line;
						line = br.readLine(); 
						continue;	
					}
					//lines that begin with ";" are comments which should be ignored
					if (line.startsWith(";")) {
						line = br.readLine(); 
						continue;	
					}
					sequence.append(line);
					line = br.readLine();
				}
				String acidSequence = sequence.toString().toUpperCase();
				nucleotideSequences.add(new DNA_Sequence(sequenceDescription, acidSequence, this));
				
				//close out our stream.  It's the courteous thing to do!
				br.close();
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			return nucleotideSequences;
		} else {
			return nucleotideSequences;
		}
	}


	
	//Resets the sequence as if it is being read for the first time
	//clears nucleotide data
	//resets index variables
	public void reset() {
		if (nucleotideSequences != null) {
			nucleotideSequences.clear();
		}
		nucleotideSequences = null;
		startIndex = 0;
		stopIndex = 0;
		nucleotideSequenceIndex = 0;
	}

}
