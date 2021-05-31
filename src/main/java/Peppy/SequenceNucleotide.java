package Peppy;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;


/**
 * A "sequence" is a DNA sequence, typically in a FASTA file
 * This farms out the digestion of the nucleotide sequences to
 * SequenceDigestionTread objects
 * 
 * Copyright 2013, Brian Risk
 * 
 * 
 * 
 * @author Brian Risk
 *
 */
public class SequenceNucleotide extends Sequence{

	private ArrayList<NucleotideSequence> nucleotideSequences = null;
	private int startIndex = 0;
	private int stopIndex = 0;
	private int nucleotideSequenceIndex = 0;
	//To avoid missing some peptides when digestion in discrete windows we
	//are allowing for some overlap
	private final int digestionFrameOverlap = 120
			; //divisible by three
	
	public SequenceNucleotide(File sequenceFile) {
		this.sequenceFile = sequenceFile;
	}
	
	public SequenceNucleotide(String sequenceFileName) {
		this(new File(sequenceFileName));
	}
	
	/**
	 * The amount of peptides produced by a sequence may be huge.
	 * therefore, this method is set up so that it should be called
	 * multiple times, each time it is called it returns another portion of the
	 * database of peptides.  When there are no more peptides from that sequence
	 * then null is returned.
	 * @param isReverse This is if we want to reverse our database to get a null database
	 * @return a sorted ArrayList of amino acid sequence fragments from the given sequence file
	 */
	public ArrayList<Peptide> extractMorePeptides(boolean isReverse) {
		//if we have yet to read in all the ATGC data, do that now
		if (nucleotideSequences == null) getNucleotideSequences();
		
		/* initialize where we will hold the peptides */
		ArrayList<Peptide> peptides = new ArrayList<Peptide>();
		
		/* we may be dealing with a set of small nucleotide sequences (in the case of RNA, cDNA)
		 * these are digested in a different way to be more efficient
		 */
		if (nucleotideSequences.get(0).getSequence().length() <= 1000000) {
			if (nucleotideSequenceIndex != -1) {
				ShortNucleotideDigestionServer snds = new ShortNucleotideDigestionServer(nucleotideSequences, isReverse);
				ProteinDigestionServer pds = new ProteinDigestionServer(snds.getResults());
				peptides = pds.getPeptides();
				/* reduce this to one peptide for a given peptide sequence */
				Hashtable<String, Peptide> uniquePeptides = new Hashtable<String, Peptide>();
				for (Peptide peptide: peptides) {
					uniquePeptides.put(peptide.getAcidSequenceString(), peptide);
				}
				peptides = new ArrayList<Peptide>(uniquePeptides.values());
				
				/* this -1 send a signal for if this function is called again that everything is already digested */
				nucleotideSequenceIndex = -1;
			} else {
				return null;
			}
		} else {
		
			//Get whatever sequence we're working on.  There will most often only be
			//one sequence in nucleotideSequences, but sometimes more
			NucleotideSequence nucleotideSequence = nucleotideSequences.get(nucleotideSequenceIndex);
			
			
			
			/* this will help us know when to clear our nucleotide data */
			boolean atLastBit = false;
			
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
							clearNucleotideSequences();
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
				if (stopIndex >= nucleotideSequence.getSequence().length()) {
					stopIndex = nucleotideSequence.getSequence().length();
					if (nucleotideSequences.size() == nucleotideSequenceIndex - 1) {
						atLastBit = true;
					}
				}
			} else {
				//if we've already extracted the region
				if (stopIndex == Properties.sequenceRegionStop) {
					//signal that we've finished
					clearNucleotideSequences();
					return null;
				} else {
					startIndex = Properties.sequenceRegionStart;
					stopIndex = Properties.sequenceRegionStop;
					atLastBit = true;
				}
			}
			
			/* Create our SequenceDigestionThread ArrayList */
			ArrayList<DigestionThread_DNA> digestors = new ArrayList<DigestionThread_DNA>();
			for (int frame = 0; frame < 3; frame++) {
				digestors.add(new DigestionThread_DNA(nucleotideSequence, frame, true, startIndex - digestionFrameOverlap, stopIndex, isReverse));
				if (!Properties.useOnlyForwardsFrames) {
					digestors.add(new DigestionThread_DNA(nucleotideSequence, frame, false, startIndex - digestionFrameOverlap, stopIndex, isReverse));
				}
			}
			
			/* create the threads and start them engines! */
			ArrayList<Thread> threads = new ArrayList<Thread>(); 
			for (DigestionThread_DNA digestor: digestors) {
				Thread thread = new Thread(digestor);
				thread.start();
				threads.add(thread);
			}
			
			/* Wait for them all to finish */
			for (Thread thread: threads) {
				try {
					thread.join();
				} catch (InterruptedException e) {
					U.p("Digestion thread interrupted!  Bad!");
					e.printStackTrace();
				}
			}
			
			/* harvest all digested peptides */
			for (DigestionThread_DNA digestor: digestors) {
				peptides.addAll(digestor.getPeptides());
			}
			
			/* free up nucleotides */
			if(atLastBit) {
				clearNucleotideSequences();
			}
		}
		
		/* sort our list of peptides by mass */
		Collections.sort(peptides);
		
		/* return */
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
	
	public ArrayList<Peptide> extractPeptidesFromRegion(int startIndex, int stopIndex, boolean isReverse) {
		
		//if we have yet to read in all the ATGC data, do that now
		if (nucleotideSequences == null) getNucleotideSequences();
		
		//Get whatever sequence we're working on.  There will most often only be
		//one sequence in nucleotideSequences, but sometimes more
		NucleotideSequence nucleotideSequence = nucleotideSequences.get(nucleotideSequenceIndex);
		
		/* bounds checking */
		if (startIndex < 0) startIndex = 0;
		if (stopIndex > nucleotideSequence.getSequence().length()) stopIndex = nucleotideSequence.getSequence().length();
		
		/* initialize where we will hold the peptides */
		ArrayList<Peptide> peptides = new ArrayList<Peptide>();
		
	
		/* Create our SequenceDigestionThread ArrayList */
		ArrayList<DigestionThread_DNA> digestors = new ArrayList<DigestionThread_DNA>();
		for (byte frame = 0; frame < 3; frame++) {
			digestors.add(new DigestionThread_DNA(nucleotideSequence, frame, true, startIndex, stopIndex, isReverse));
			if (!Properties.useOnlyForwardsFrames) {
				digestors.add(new DigestionThread_DNA(nucleotideSequence, frame, false, startIndex, stopIndex, isReverse));
			}
		}
		
		/* create the threads and start them engines! */
		ArrayList<Thread> threads = new ArrayList<Thread>(); 
		for (DigestionThread_DNA digestor: digestors) {
			Thread thread = new Thread(digestor);
			thread.start();
			threads.add(thread);
		}
		
		/* Wait for them all to finish */
		for (Thread thread: threads) {
			try {
				thread.join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		
		/* harvest all digested peptides */
		for (DigestionThread_DNA digestor: digestors) {
			peptides.addAll(digestor.getPeptides());
		}
		
		/* return */
		return peptides;
	}
	
	
	/**
	 * in one FASTA file there may be many sequences.
	 * this method returns a ArrayList containing all of the sequences
	 * @return
	 */
	public ArrayList<NucleotideSequence> getNucleotideSequences() {
		if (nucleotideSequences == null) {
			nucleotideSequences = new ArrayList<NucleotideSequence>();
			try {
				BufferedReader br = new BufferedReader(new FileReader(sequenceFile));
				String line = br.readLine();

				String sequenceDescription = null;
				StringBuffer sequence = new StringBuffer();
				while (line != null) {
					
					/* lines that begin with ">" are comments which name the sequence */
					if (line.startsWith(">")) {
						
						/* if sequenceDescription has not been defined, that means it is the first in the file */
						if (sequenceDescription != null) {
							nucleotideSequences.add(new NucleotideSequence(sequenceDescription, sequence.toString().toUpperCase(), this));
							sequence = new StringBuffer();
						}
						sequenceDescription = line;
						line = br.readLine(); 
						continue;	
					}
					
					/* lines that begin with ";" are comments which should be ignored */
					if (line.startsWith(";")) {
						line = br.readLine(); 
						continue;	
					}
					sequence.append(line);
					line = br.readLine();
				}
				String acidSequence = sequence.toString().toUpperCase();
//				String acidSequence = sequence.reverse().toString().toUpperCase();
				nucleotideSequences.add(new NucleotideSequence(sequenceDescription, acidSequence, this));
				
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
		clearNucleotideSequences();
		nucleotideSequences = null;
		startIndex = 0;
		stopIndex = 0;
		nucleotideSequenceIndex = 0;
	}

	private void clearNucleotideSequences() {
		if (nucleotideSequences != null) {
			for (NucleotideSequence seq: nucleotideSequences) {
				seq.clearSequenceData();
			}
			nucleotideSequences.clear();
		}
//		System.gc();
	}

}
