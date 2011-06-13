package Peppy;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;

import Utilities.U;

/**
 * The contract here is that all lists of peptides must be returned
 * sorted by mass, from least to greatest.
 * @author Brian Risk
 *
 */
public class Sequence_Protein extends Sequence {
	
	BufferedReader reader = null;
	ArrayList<Protein> proteins = null;
	
	public Sequence_Protein(File proteinFile) {
		this.sequenceFile = proteinFile;
		try {
			reader = new BufferedReader(new FileReader(proteinFile));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	public void reset() {
		proteins = null;
		try {
			reader = new BufferedReader(new FileReader(sequenceFile));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	public  ArrayList<Peptide> extractAllPeptides(boolean isReverse) {
		getProteinsFromDatabase(isReverse, false);
		return getPeptidesFromListOfProteins(proteins);	
	}
	
	public  ArrayList<Peptide> extractMorePeptides(boolean isReverse) {
		getProteinsFromDatabase(isReverse, true);
		return getPeptidesFromListOfProteins(proteins);	
	}
	


	/**
	 * This method, right here, it's for reverse database searches for things like
	 * False Positive Rates and the like.  So, only used for stats purposes, not
	 * "real" searches
	 * @param proteinFile
	 * @return
	 */
	public  ArrayList<Peptide> getAllPeptidesFromReverseDatabase(File proteinFile) {
		getProteinsFromDatabase(true, false);
		return getPeptidesFromListOfProteins(proteins);	
	}
	
	/**
	 * Depending on the file suffix of the protein file it chooses how to extract the proteins
	 * @param proteinFile a FASTA or DAT formatted file
	 * @param isReverse
	 * @return
	 */
	private ArrayList<Protein> getProteinsFromDatabase(boolean isReverse, boolean limitAmount) {
		if (sequenceFile.getName().toLowerCase().endsWith("fa")) return getProteinsFromFASTA(isReverse, limitAmount);
		if (sequenceFile.getName().toLowerCase().endsWith("fsa")) return getProteinsFromFASTA(isReverse, limitAmount);
		if (sequenceFile.getName().toLowerCase().endsWith("fasta")) return getProteinsFromFASTA(isReverse, limitAmount);
		if (sequenceFile.getName().toLowerCase().endsWith("txt")) return getProteinsFromFASTA(isReverse, limitAmount);
		if (sequenceFile.getName().toLowerCase().endsWith("dat")) return getProteinsFromUniprotDAT(isReverse, limitAmount);
		return null;
	}
	
	public static ArrayList<Peptide> getPeptidesFromListOfProteins(ArrayList<Protein> proteins) {
		if (proteins.size() == 0) return null;
		ArrayList<Peptide> peptides = new ArrayList<Peptide>();
		for (Protein protein: proteins) {
			peptides.addAll(protein.digest());
		}
		Collections.sort(peptides);
		return peptides;
	}
	
	public static void main(String [] args) {
		String test = "piowenrweoiren*sflkwnreoin*powiernewoinr*oin";
		String [] chunks = test.split("\\*");
		for (int i = 0; i < chunks.length; i++) {
			U.p(chunks[i]);
		}
	}
	
	private ArrayList<Protein> getProteinsFromFASTA( boolean isReverse, boolean limitAmount) {
		proteins = new ArrayList<Protein>();
		try {
			String line = reader.readLine();
			
			/*since proteins occur on multiple lines, this buffer is where they are all put together */
			StringBuffer buffy = new StringBuffer();

			/* sometimes there are "*" in protein strings to represent chunks.  
			 * this is where we store those chunks
			 */
			String [] proteinChunks;
			
			/* the accession number */
			String proteinName = "";
			while (line != null) {
				
				/* this symbol means we've reached the beginning of a new protein and
				the one we've been working on has ended */
				if (line.startsWith(">")) {
					
					/* make a new protein if we've been building one */
					if (buffy.length() > 0) {
						
						/* reverse the string if we are making a null database */
						if (isReverse) buffy.reverse();

						/* again, this addresses the "*" in some FASTA files */
						proteinChunks = buffy.toString().split("\\*");
						
						/* go through each of our chunks and add them as proteins */
						for (int i = 0; i < proteinChunks.length; i++) {
							proteins.add(new Protein(proteinName, buffy.toString()));
						}
						
					}

					/* this is the name of the next protein */
					proteinName = line.substring(1).trim();
					
					/* try to get the accession number from UniProt databases */
					if (proteinName.startsWith("sp|") || proteinName.startsWith("tr|")) {
						String [] proteinWords = proteinName.split("\\|");
						proteinName = proteinWords[1];
					}
					
					buffy = new StringBuffer(); 
				} else {
					buffy.append(line);
				}
				if (limitAmount) {
					if (proteins.size() > Properties.maxNumberOfProteinsToLoadAtOnce) return proteins;
				}
				line = reader.readLine();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return proteins;
	}
	
	/**
	 * Uniprot has its own format!  Yes!
	 * @param proteinFile
	 * @return
	 */
	private ArrayList<Protein> getProteinsFromUniprotDAT(boolean reverse, boolean limitAmount) {
		proteins = new ArrayList<Protein>();
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
					if (reverse) buffy.reverse();
					proteins.add(new Protein(proteinName, buffy.toString()));
				}
				if (limitAmount) {
					if (proteins.size() > Properties.maxNumberOfProteinsToLoadAtOnce) return proteins;
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
