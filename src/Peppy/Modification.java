package Peppy;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

public class Modification implements Comparable<Modification>{

	int accessionNumber;
	String PSI_MSname;
	String interimName;
	String description;
	double monoMass;
	double averageMass;
	String compositioin;
	
	/**
	 * builds a ProteinModification from string from the ptm file
	 * @param line
	 */
	public Modification(String line) {
		String [] chunks = line.split("\t");
		accessionNumber = Integer.parseInt(chunks[0]);
		PSI_MSname = chunks[1];
		interimName = chunks[2];
		description = chunks[3];
		monoMass = Double.parseDouble(chunks[4]);
		averageMass = Double.parseDouble(chunks[5]);
		compositioin = chunks[6];
	}
	
	public static ArrayList<Modification> getProteinModificationsFromFile(File file) {
		ArrayList<Modification> out = new ArrayList<Modification>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			String line = br.readLine();
			while (line != null) {
				out.add(new Modification(line));
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

	public int compareTo(Modification other) {
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
