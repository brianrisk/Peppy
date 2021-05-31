package Tools.VCFConverter;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class VCFFile {
	ArrayList<VCFEntry> entries = new ArrayList<VCFEntry>();
	ArrayList<VCFEntry> invalidEntries = new ArrayList<VCFEntry>();
	
	public VCFFile (String file) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			String line = br.readLine();
			line = br.readLine();
			while (line != null) {
				if (line.startsWith("#")) {
					line = br.readLine();
					continue;
				}
				VCFEntry ve = new VCFEntry(line);
				if (ve.isValid()){
					entries.add(ve);
				}
				line = br.readLine();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(1);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public ArrayList<VCFEntry> getEntries() {
		return entries;
	}

	public ArrayList<VCFEntry> getInvalidEntries() {
		return invalidEntries;
	}

}
