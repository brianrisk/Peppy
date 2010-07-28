package Peppy;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collections;
import java.util.ArrayList;

import Utilities.U;

/**
 * A "sequence" is a DNA sequence, typically in a FASTA file
 * This farms out the digestion of the nucleotide sequences to
 * SequenceDigestionTread objects
 * @author Brian Risk
 *
 */
public class Sequence {
	private File sequenceFile;
	private int id;
	private ArrayList<NucleotideSequence> nucleotideSequences = null;
	private int startIndex = 0;
	private int stopIndex = 0;
	private int digestionWindowSize = 25000000;
	private int nucleotideSequenceIndex = 0;
	
	public Sequence(File sequenceFile) {
		this.sequenceFile = sequenceFile;
	}
	
	public Sequence(String sequenceFileName) {
		this(new File(sequenceFileName));
	}
	
	/**
	 * The amount of peptides produced by a sequence may be huge.
	 * therefore, this method is set up so that it should be called
	 * multiple times, each time it is called it returns another portion of the
	 * database of peptides.
	 * @return a sorted ArrayList of amino acid sequence fragments from the given sequence file
	 */
	public ArrayList<Peptide> extractPeptides() {
		if (nucleotideSequences == null) getNucleotideSequences();
		NucleotideSequence nucleotideSequence = nucleotideSequences.get(nucleotideSequenceIndex);
		
		//if this is the first time we've tried to extract peptides
		if (stopIndex == 0) {
			stopIndex = digestionWindowSize;
		} else {
		
			//if we've reached the end of the nucleotide sequence
			if (stopIndex == nucleotideSequence.getSequence().length()) {
				nucleotideSequenceIndex++;
				if (nucleotideSequenceIndex >= nucleotideSequences.size()) return null;
				nucleotideSequence = nucleotideSequences.get(nucleotideSequenceIndex);
				startIndex = 0;
				stopIndex = digestionWindowSize;
			} else {
				startIndex += digestionWindowSize;
				stopIndex += digestionWindowSize;
			}
		}
		//overshoot correction
		if (stopIndex > nucleotideSequence.getSequence().length()) {
			stopIndex = nucleotideSequence.getSequence().length();
		}
		
		
		ArrayList<Peptide> peptides = new ArrayList<Peptide>();
		
		//Create our SequenceDigestionThread ArrayList
		ArrayList<SequenceDigestionThread> digestors = new ArrayList<SequenceDigestionThread>();
		for (byte frame = 0; frame < 3; frame++) {
			digestors.add(new SequenceDigestionThread(nucleotideSequence, frame, true, startIndex, stopIndex));
			digestors.add(new SequenceDigestionThread(nucleotideSequence, frame, false, startIndex, stopIndex));
		}
		
		ArrayList<Thread> threads = new ArrayList<Thread>();
		for (int digestorIndex = 0; digestorIndex < digestors.size(); digestorIndex++) {
			Thread thread = new Thread(digestors.get(digestorIndex));
			thread.start();
			threads.add(thread);
		}
		
		//Wait for them all to finish
		for (int threadIndex = 0; threadIndex < threads.size(); threadIndex++) {
			Thread thread = threads.get(threadIndex);
			try {
				thread.join();
			} catch (InterruptedException e) {
				U.p("Digestion thread interrupted!  Bad!");
				e.printStackTrace();
			}
		}
		
		//harvest all digested peptides
		for (int digestorIndex = 0; digestorIndex < digestors.size(); digestorIndex++) {
			SequenceDigestionThread digestor = digestors.get(digestorIndex);
			peptides.addAll(digestor.getPeptides());
		}
			
		Collections.sort(peptides);
		return peptides;
	}
	

	public static ArrayList<Sequence> loadSequences(File folder) { 
		ArrayList<Sequence> sequences = new ArrayList<Sequence>();
		if (folder.isFile()) {
			sequences.add(new Sequence(folder));
		} else {
			File [] files = folder.listFiles();
			for (int i = 0; i < files.length; i++) {
				if (files[i].isHidden()) continue;
				String fileName = files[i].getName().toLowerCase();
				if (fileName.endsWith(".fasta") || fileName.endsWith(".fa") || fileName.endsWith(".fsa") || fileName.endsWith(".txt")) {
					sequences.add(new Sequence(files[i]));
				}
			}
		}
		for (int i = 0; i < sequences.size(); i++) {
			sequences.get(i).setId(i);
		}
		return sequences;
	}
	
	public static ArrayList<Sequence> loadSequences(String folderName) { 
		return loadSequences(new File(folderName));
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
				//lines beginning with ">" name the sequence
				String sequenceDescription = null;
				StringBuffer sequence = new StringBuffer();
				while (line != null) {
					//lines that begin with ">" are comments which name the sequence
					if (line.startsWith(">")) {
						//if sequenceDescription has not been defined, that means it is the first in the file
						if (sequenceDescription != null) {
							nucleotideSequences.add(new NucleotideSequence(sequenceDescription, new String(sequence), this));
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
					line = br.readLine().toUpperCase();
				}
				nucleotideSequences.add(new NucleotideSequence(sequenceDescription, new String(sequence), this));
				
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

	public File getSequenceFile() {
		return sequenceFile;
	}

	public int getId() {
		return id;
	}

	public void setId(int id) {
		this.id = id;
	}
	
	//TODO
	/**
	 * Since there can be many sequences inside of one sequence file this
	 * method is not the best.  It needs some sort of improvement.  In fact,
	 * the whole basic data structure for this kind of thing may need to be
	 * overhauled.
	 * 
	 * Oh well.  In the mean time this returns the length of the first
	 * nucleotide sequence.
	 * 
	 */
	public int getSequenceLength() {
		if (nucleotideSequences == null) {
			return 0;
		} else {
			if (nucleotideSequences.size() == 0) {
				return 0;
			} else {
				NucleotideSequence ns = nucleotideSequences.get(0);
				return ns.getSequence().length();
			}
		}
	}

}
