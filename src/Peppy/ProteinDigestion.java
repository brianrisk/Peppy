package Peppy;
import java.io.*;
import java.util.ArrayList;
import java.util.Collections;

import Utilities.U;


public class ProteinDigestion {
	
	public static ArrayList<Peptide> getPeptidesFromDatabase(File proteinFile) {
		return getPeptidesFromDatabase(proteinFile, false);
	}
	
	public static ArrayList<Peptide> getPeptidesFromDatabase(File proteinFile, boolean reverse) {
		if (proteinFile.getName().toLowerCase().endsWith("fa")) return getPeptidesFromFASTA(proteinFile, reverse);
		if (proteinFile.getName().toLowerCase().endsWith("fsa")) return getPeptidesFromFASTA(proteinFile, reverse);
		if (proteinFile.getName().toLowerCase().endsWith("fasta")) return getPeptidesFromFASTA(proteinFile, reverse);
		if (proteinFile.getName().toLowerCase().endsWith("dat")) return getPeptidesFromUniprotDAT(proteinFile, reverse);
		return null;
	}
	
	private static ArrayList<Peptide> getPeptidesFromFASTA(File proteinFile, boolean reverse) {
		ArrayList<Peptide> out = new ArrayList<Peptide>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(proteinFile));
			String line = br.readLine();
			StringBuffer buffy = new StringBuffer();
			String proteinName = "";
			int proteinIndex = 0;
			Protein protein;
			while (line != null) {
				if (line.startsWith(">")) {
					if (reverse) buffy.reverse();
					protein = new Protein(proteinName, buffy.toString());
					out.addAll(protein.digest());
					proteinName = line.substring(1).trim();
					//try to get the accession number from UniProt databases
					if (proteinName.startsWith("sp|") || proteinName.startsWith("tr|")) {
						String [] proteinWords = proteinName.split("|");
						proteinName = proteinWords[1];
					}
					
					buffy = new StringBuffer(); 
					proteinIndex++;
				} else {
					buffy.append(line);
				}
				line = br.readLine();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		Collections.sort(out);
		return out;
	}
	
	/**
	 * Uniprot has its own format!  Yes!
	 * @param proteinFile
	 * @return
	 */
	private static ArrayList<Peptide> getPeptidesFromUniprotDAT(File proteinFile, boolean reverse) {
		ArrayList<Peptide> out = new ArrayList<Peptide>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(proteinFile));
			String line = br.readLine();
			StringBuffer buffy = new StringBuffer();
			String proteinName = "";
			boolean inSequence = false;
			Protein protein;
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
					protein = new Protein(proteinName, buffy.toString());
					out.addAll(protein.digest());
				}
				//read a new line
				line = br.readLine();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		Collections.sort(out);
		return out;
	}
	
	//TODO
	/**
	 * DELETE LATER
	 * this is a quick converter from DAT to FASTA
	 * @param proteinFile
	 * @return
	 */
	public static void convertDATtoFASTA(File proteinFile) {
		U.p("loading protein file: " + proteinFile.getName());
		try {
			BufferedReader br = new BufferedReader(new FileReader(proteinFile));
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File("uniprot_human.fasta"))));
			String line = br.readLine();
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
					pw.println("> " + proteinName);
					String protein = buffy.toString();
					for (int i = 0; i < protein.length(); i++) {
						pw.print(protein.charAt(i));
						if (((i + 1) % 60 == 0) && i > 1) {
							pw.print("\r");
						}
					}
					pw.print("\r");
					
				}
				//read a new line
				line = br.readLine();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	
	

}
