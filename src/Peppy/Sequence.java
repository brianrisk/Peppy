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
	
	public Sequence(File sequenceFile) {
		this.sequenceFile = sequenceFile;
	}
	
	public Sequence(String sequenceFileName) {
		this(new File(sequenceFileName));
	}
	
	/**
	 * 
	 * @return an sorted ArrayList of amino acid sequence fragments from the given sequence file
	 */
	public ArrayList<Peptide> extractPeptides() {
		U.p("digesting sequences...");
		ArrayList<NucleotideSequence> nucleotideSequences = getNucleotideSequences();
		ArrayList<Peptide> peptides = new ArrayList<Peptide>();
		for (int sequenceIndex = 0; sequenceIndex < nucleotideSequences.size(); sequenceIndex++) {
			NucleotideSequence nucleotideSequence = nucleotideSequences.get(sequenceIndex);
			
			//Create our SequenceDigestionThread ArrayList
			ArrayList<SequenceDigestionThread> digestors = new ArrayList<SequenceDigestionThread>();
			for (byte frame = 0; frame < 3; frame++) {
				digestors.add(new SequenceDigestionThread(nucleotideSequence, frame, true, false));
				digestors.add(new SequenceDigestionThread(nucleotideSequence, frame, false, false));
				//digestors.add(new SequenceDigestionThread(nucleotideSequence, frame, true, true));
				//digestors.add(new SequenceDigestionThread(nucleotideSequence, frame, false, true));
			}
			
			//Putting the thread ArrayList here gives us a max of 12 threads.
			//put it outside the other "for" if you have multiple sequences
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
			
		}
		Collections.sort(peptides);
		U.p("done digesting.");
		return peptides;
	}
	
	/**
	 * /**
	 * This method, though it uses a threads, is not really for multi threading.
	 * It is for when we want to be memory-conscious and only digest a sequence one way
	 * (reading frame / missed cleavage) at a time.
	 * @param nucleotideSequences
	 * @param frame
	 * @param forwards
	 * @param missedCleavage
	 * @return an sorted ArrayList of amino acid sequence fragments from the given nucleotide sequence
	 */
	public ArrayList<Peptide> extractPeptides(ArrayList<NucleotideSequence> nucleotideSequences, byte frame, boolean forwards, boolean missedCleavage) {
		ArrayList<Peptide> peptides = null;
		for (int sequenceIndex = 0; sequenceIndex < nucleotideSequences.size(); sequenceIndex++) {
			NucleotideSequence nucleotideSequence = nucleotideSequences.get(sequenceIndex);
			SequenceDigestionThread digestor = new SequenceDigestionThread(nucleotideSequence, frame, forwards, missedCleavage);
			Thread thread = new Thread(digestor);
			thread.start();
			try {
				thread.join();
			} catch (InterruptedException e) {
				U.p("Digestion thread interrupted!  Bad!");
				e.printStackTrace();
			}
			peptides = digestor.getPeptides();	
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
				if (fileName.endsWith(".fasta") || fileName.endsWith(".fa")) {
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
		ArrayList<NucleotideSequence> out = new ArrayList<NucleotideSequence>();
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
						out.add(new NucleotideSequence(sequenceDescription, new String(sequence), this));
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
			out.add(new NucleotideSequence(sequenceDescription, new String(sequence), this));
			
			//close out our stream.  It's the courteous thing to do!
			br.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return out;
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

}
