package Peppy;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class ModificationEntry implements Comparable<ModificationEntry>{

	int accessionNumber;
	String PSI_MSname;
	String interimName;
	String description;
	double monoMass;
	double averageMass;
	String compositioin;
	
	//no modification
	public ModificationEntry() {
		accessionNumber = -1;
		PSI_MSname = "none";
		interimName = "none";
		description = "none";
		monoMass = 0;
		averageMass = 0;
		compositioin = "";
	}
	
	/**
	 * builds a ProteinModification from string from the ptm file
	 * @param line
	 */
	public ModificationEntry(String line) {
		String [] chunks = line.split("\t");
		accessionNumber = Integer.parseInt(chunks[0]);
		PSI_MSname = chunks[1];
		interimName = chunks[2];
		description = chunks[3];
		monoMass = Double.parseDouble(chunks[4]);
		averageMass = Double.parseDouble(chunks[5]);
		compositioin = chunks[6];
	}
	
	public static ArrayList<ModificationEntry> getProteinModificationsFromFile(File file) {
		ArrayList<ModificationEntry> out = new ArrayList<ModificationEntry>();
		//add the null modification
		out.add(new ModificationEntry());
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			String line = br.readLine();
			while (line != null) {
				out.add(new ModificationEntry(line));
				line = br.readLine();
			}
			br.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		Collections.sort(out);
		return out;
	}

	public int compareTo(ModificationEntry other) {
		if (getMonoMass() < other.getMonoMass()) return -1;
		if (getMonoMass() > other.getMonoMass()) return  1;
		return 0;
	}
	
	public String toString() {
		return description;
	}
	
	public int getAccessionNumber() {
		return accessionNumber;
	}

	public String getPSI_MSname() {
		return PSI_MSname;
	}

	public String getInterimName() {
		return interimName;
	}

	public String getDescription() {
		return description;
	}

	public double getMonoMass() {
		return monoMass;
	}

	public double getAverageMass() {
		return averageMass;
	}

	public String getCompositioin() {
		return compositioin;
	}
	

}
