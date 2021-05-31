package Peppy;

import java.io.*;
import java.util.ArrayList;

/**
 * The contract here is that all lists of peptides must be returned
 * sorted by mass, from least to greatest.
 * 
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class SequenceProtein extends Sequence {
	
	private BufferedReader reader = null;
	
	public SequenceProtein(File proteinFile) {
		this.sequenceFile = proteinFile;
		try {
			reader = new BufferedReader(new FileReader(proteinFile));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	public void reset() {
		try {
			reader = new BufferedReader(new FileReader(sequenceFile));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	public  ArrayList<Peptide> extractAllPeptides(boolean isReverse) {
		ArrayList<Protein> proteins = getProteinsFromDatabase(isReverse, false);
		return getPeptidesFromListOfProteins(proteins);
	}
	
	public  ArrayList<Peptide> extractMorePeptides(boolean isReverse) {
		ArrayList<Protein> proteins = getProteinsFromDatabase(isReverse, true);
		ArrayList<Peptide> peptides = getPeptidesFromListOfProteins(proteins);
		if (peptides.size() == 0) return null; //TODO fix.  there is the (small) possibility that this may return null and there are still proteins left to digest.
		return peptides;
	}
	


	/**
	 * This method, right here, it's for reverse database searches for things like
	 * False Positive Rates and the like.  So, only used for stats purposes, not
	 * "real" searches
	 * @param proteinFile
	 * @return
	 */
	public  ArrayList<Peptide> getAllPeptidesFromReverseDatabase(File proteinFile) {
		ArrayList<Protein> proteins = getProteinsFromDatabase(true, false);
		return getPeptidesFromListOfProteins(proteins);
	}
	
	/**
	 * Depending on the file suffix of the protein file it chooses how to extract the proteins
	 * @param isReverse
	 * @return
	 */
	public ArrayList<Protein> getProteinsFromDatabase(boolean isReverse, boolean limitAmount) {
		if (sequenceFile.getName().toLowerCase().endsWith("fa")) return getProteinsFromFASTA(isReverse, limitAmount);
		if (sequenceFile.getName().toLowerCase().endsWith("fsa")) return getProteinsFromFASTA(isReverse, limitAmount);
		if (sequenceFile.getName().toLowerCase().endsWith("fasta")) return getProteinsFromFASTA(isReverse, limitAmount);
		if (sequenceFile.getName().toLowerCase().endsWith("txt")) return getProteinsFromFASTA(isReverse, limitAmount);
		if (sequenceFile.getName().toLowerCase().endsWith("dat")) return getProteinsFromUniprotDAT(isReverse, limitAmount);
		return null;
	}
	
	public static ArrayList<Peptide> getPeptidesFromListOfProteins(ArrayList<Protein> proteins) {
		ArrayList<Peptide> peptides = new ArrayList<>();
		if (proteins.size() != 0) {
			ProteinDigestionServer pds = new ProteinDigestionServer(proteins);
			peptides = pds.getPeptides();
		}
		return peptides;
	}
	
	private ArrayList<Protein> getProteinsFromFASTA( boolean isReverse, boolean limitAmount) {
		ArrayList<Protein> proteins = new ArrayList<>();
		int combinedLength = 0;
		try {
			String line = reader.readLine();
			
			/*since proteins occur on multiple lines, this buffer is where they are all put together */
			StringBuffer buffy = new StringBuffer();
			
			/* the accession number */
			String proteinName = "";
			while (line != null) {
				
				/* this symbol means we've reached the beginning of a new protein and
				the one we've been working on has ended */
				if (line.startsWith(">")) {
				
					/* keep track of length */
					combinedLength += buffy.length();
					
					/* make a new protein if we've been building one */
					addProtein(buffy, proteinName, proteins, isReverse);	

					/* this is the name of the next protein */
					proteinName = line.substring(1).trim();
					
					/* try to get the accession number from UniProt databases */
					if (proteinName.startsWith("sp|") || proteinName.startsWith("tr|")) {
						String [] proteinWords = proteinName.split("\\|");
						proteinName = proteinWords[1];
					}
					
					/* initialize our string buffer */
					buffy = new StringBuffer(); 
				} else {
					buffy.append(line);
				}
				if (limitAmount) {
					if (combinedLength > Properties.maxCombinedProteinLength) return proteins;
				}
				line = reader.readLine();
			}
			
			/* Add the final line */
			addProtein(buffy, proteinName, proteins, isReverse);
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		return proteins;
	}
	
	/**
	 * Takes a string buffer and adds it to our list of proteins with the appropriate name
	 */
	private void addProtein(StringBuffer buffy, String proteinName, ArrayList<Protein> proteins, boolean isReverse) {
		if (buffy.length() > 0) {
			
			/* sometimes there are "*" in protein strings to represent chunks.  
			 * this is where we store those chunks
			 */
			String [] proteinChunks;
			
			/* reverse the string if we are making a null database */
			if (isReverse) buffy.reverse();
	
			/* this addresses the "*" in some FASTA files */
			proteinChunks = buffy.toString().toUpperCase().split("\\*");
			
			/* go through each of our chunks and add them as proteins */
			for (String proteinChunk : proteinChunks) {
				if (proteinChunk.length() < Properties.minPeptideLength) continue;
				proteins.add(new Protein(proteinName, proteinChunk, isReverse));
			}
		}
	}
	
	/**
	 * Uniprot has its own format!  Yes!
	 */
	private ArrayList<Protein> getProteinsFromUniprotDAT(boolean isDecoy, boolean limitAmount) {
		ArrayList<Protein> proteins = new ArrayList<>();
		try {
			String line;
			line = reader.readLine();
			StringBuffer buffy = new StringBuffer();
			String proteinName = "";
			boolean inSequence = false;
			while (line != null) {				
				if (line.startsWith("ID")) {
					proteinName = "NOT DEFINED";
					buffy = new StringBuffer(); 
				}
				//DR   UCSC; uc002hjd.1; human.
				if (line.startsWith("AC")) {
					if (proteinName.equals("NOT DEFINED")) {
						String [] chunks = line.split(";");
						proteinName = chunks[0].trim();
						proteinName = proteinName.substring(5);
					}
				}
				if (line.startsWith("SQ")) {
					inSequence = true;
				}
				//starts with five spaces
				if (line.startsWith("     ")) {
					if (inSequence) {
						char theChar;
						for (int i = 0; i < line.length(); i++) {
							theChar = line.charAt(i);
							if (theChar != ' ') buffy.append(theChar);
						}
					}
				}
				if (line.startsWith("//")) {
					inSequence = false;
					if (isDecoy) buffy.reverse();
					proteins.add(new Protein(proteinName, buffy.toString(), isDecoy));
				}
				if (limitAmount) {
					if (proteins.size() > 5000) return proteins;
				}
				//read a new line
				line = reader.readLine();
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return proteins;
	}
	
	
	
	

}
