package Peppy;
import java.io.*;
import java.util.ArrayList;
import java.util.Collections;

import Utilities.U;


public class ProteinDigestion {
	
	public static ArrayList<Peptide> getPeptidesFromDatabase(File proteinFile) {
		if (proteinFile.getName().toLowerCase().endsWith("fa")) return getPeptidesFromFASTA(proteinFile);
		if (proteinFile.getName().toLowerCase().endsWith("fsa")) return getPeptidesFromFASTA(proteinFile);
		if (proteinFile.getName().toLowerCase().endsWith("fasta")) return getPeptidesFromFASTA(proteinFile);
		if (proteinFile.getName().toLowerCase().endsWith("dat")) return getPeptidesFromUniprotDAT(proteinFile);
		return null;
	}
	
	public static ArrayList<Peptide> getPeptidesFromFASTA(File proteinFile) {
		ArrayList<Peptide> out = new ArrayList<Peptide>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(proteinFile));
			String line = br.readLine();
			StringBuffer buffy = new StringBuffer();
			String proteinName = "";
			while (line != null) {
				if (line.startsWith(">")) {
					out.addAll(getPeptidesFromProteinString(buffy.toString(), proteinName));
					proteinName = line;
					buffy = new StringBuffer(); 
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
	public static ArrayList<Peptide> getPeptidesFromUniprotDAT(File proteinFile) {
		ArrayList<Peptide> out = new ArrayList<Peptide>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(proteinFile));
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
					out.addAll(getPeptidesFromProteinString(buffy.toString(), proteinName));
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
	
	public static ArrayList<Peptide> getReversePeptidesFromFASTA(File proteinFile) {
		ArrayList<Peptide> out = new ArrayList<Peptide>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(proteinFile));
			String line = br.readLine();
			StringBuffer buffy = new StringBuffer();
			String proteinName = "";
			while (line != null) {
				if (line.startsWith(">")) {
					String forwards = buffy.toString();
					StringBuffer reverseBuffer = new StringBuffer();
					for (int i = forwards.length() - 1; i >=0; i--) {
						reverseBuffer.append(forwards.charAt(i));
					}
					out.addAll(getPeptidesFromProteinString(reverseBuffer.toString(), proteinName));
					buffy = new StringBuffer(); 
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
	
	public static ArrayList<Peptide> getPeptidesFromProteinString(String proteinString, String proteinName) {
		ArrayList<Peptide> out = new ArrayList<Peptide>();
		if (proteinString.length() == 0) return out;
		ArrayList<Peptide> fragments = new ArrayList<Peptide>();
		StringBuffer buffy = new StringBuffer();
		boolean cleavage;
		//starting with the second amino acid
		for (int i = 1; i < proteinString.length(); i++) {
			buffy.append(proteinString.charAt(i - 1));
			cleavage = false;
			if (proteinString.charAt(i - 1) == 'K') cleavage = true;
			if (proteinString.charAt(i - 1) == 'R') cleavage = true;
			if (proteinString.charAt(i) == 'P') cleavage = false;
			if (proteinString.charAt(i) == 'X') cleavage = true;
			if (cleavage) {
				Peptide peptide = new Peptide(buffy.toString(), proteinName);
				fragments.add(peptide);
				buffy = new StringBuffer();
			}
			if (proteinString.charAt(i) == 'X') {
				while (proteinString.charAt(i) == 'X') {
					i++;
					if (i ==  proteinString.length()) break;
				}
				i++;
			}
		}
		
		//get in the last peptide
		buffy.append(proteinString.charAt(proteinString.length() - 1));
		fragments.add(new Peptide(buffy.toString(), proteinName));
		
		//add big enough fragments to out
		for (int i = 0; i < fragments.size(); i++) {
			Peptide peptide = fragments.get(i);
			if (peptide.getMass() >= Properties.peptideMassThreshold) out.add(peptide);
		}
		
		//getting all missed cleavages
		for (int numberOfMissedCleavages = 1; numberOfMissedCleavages <= Properties.numberOfMissedCleavages; numberOfMissedCleavages++){
			for (int i = 0; i < fragments.size() - numberOfMissedCleavages; i++) {
				StringBuffer peptideString = new StringBuffer(fragments.get(i).getAcidSequence());
				for (int j = 1; j <= numberOfMissedCleavages; j++) {
					peptideString.append(fragments.get(i + j).getAcidSequence());
				}
				Peptide peptide = new Peptide(peptideString.toString(),  proteinName);
				if (peptide.getMass() >= Properties.peptideMassThreshold) out.add(peptide);
			}
		}
		return out;
	}

}
