package Peppy;
import java.io.*;
import java.util.ArrayList;
import java.util.Collections;

import Utilities.U;


public class ProteinDigestion {
	
	public static ArrayList<Peptide> getPeptidesFromProteinFile(File proteinFile) {
		U.p("loading protein file: " + proteinFile.getName());
		ArrayList<Peptide> out = new ArrayList<Peptide>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(proteinFile));
			String line = br.readLine();
			StringBuffer buffy = new StringBuffer();
			while (line != null) {
				if (line.startsWith(">")) {
					out.addAll(getPeptidesFromProteinString(buffy.toString(), 0, true, (byte) 0, null));
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
	
	public static ArrayList<Peptide> getReversePeptidesFromProteinFile(File proteinFile) {
		U.p("loading reverse protein file: " + proteinFile.getName());
		ArrayList<Peptide> out = new ArrayList<Peptide>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(proteinFile));
			String line = br.readLine();
			StringBuffer buffy = new StringBuffer();
			while (line != null) {
				if (line.startsWith(">")) {
					String forwards = buffy.toString();
					StringBuffer reverseBuffer = new StringBuffer();
					for (int i = forwards.length() - 1; i >=0; i--) {
						reverseBuffer.append(forwards.charAt(i));
					}
					out.addAll(getPeptidesFromProteinString(reverseBuffer.toString(), 0, true, (byte) 0, null));
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
	
	public static ArrayList<Peptide> getPeptidesFromProteinString(String proteinString, int startIndex, boolean forwards, byte frame, Sequence parentSequence) {
		ArrayList<Peptide> out = new ArrayList<Peptide>();
		if (proteinString.length() == 0) return out;
		ArrayList<Peptide> fragments = new ArrayList<Peptide>();
		StringBuffer buffy = new StringBuffer();
		boolean cleavage;
		int locus;
		//starting with the second amino acid
		for (int i = 1; i < proteinString.length(); i++) {
			buffy.append(proteinString.charAt(i - 1));
			cleavage = false;
			if (proteinString.charAt(i - 1) == 'K') cleavage = true;
			if (proteinString.charAt(i - 1) == 'R') cleavage = true;
			if (proteinString.charAt(i) == 'P') cleavage = false;
			if (proteinString.charAt(i) == 'X') cleavage = true;
			if (cleavage) {
				if (forwards) locus = startIndex + ((i -1) * 3);
				else locus = startIndex - ((i -1) * 3);
				Peptide peptide = new Peptide(buffy.toString(), locus, forwards, frame, parentSequence);
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
		if (forwards) locus = startIndex + ((proteinString.length() -1) * 3);
		else locus = startIndex - ((proteinString.length() -1) * 3);
		buffy.append(proteinString.charAt(proteinString.length() - 1));
		fragments.add(new Peptide(buffy.toString(), locus, forwards, frame, parentSequence));
		
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
				Peptide peptide = new Peptide(peptideString.toString(),  fragments.get(i).getIndex(), forwards, frame, parentSequence);
				if (peptide.getMass() >= Properties.peptideMassThreshold) out.add(peptide);
			}
		}
		
		//two missed cleavages
//		for (int i = 0; i < fragments.size() - 2; i++) {
//			String peptideString = fragments.get(i).getSequence() + fragments.get(i + 1).getSequence() + fragments.get(i + 2).getSequence();
//			Peptide peptide = new Peptide(peptideString,  fragments.get(i).getIndex(), forwards, frame, parentSequence);
//			if (peptide.getMass() >= Properties.peptideMassThreshold) out.add(peptide);
//		}
		return out;
	}

}
