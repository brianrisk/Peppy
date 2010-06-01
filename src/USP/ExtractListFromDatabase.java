package USP;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import Utilities.U;

/**
 * This is the first step in the USP evaluation.  We are first
 * taking the list of proteins from United Proteomics Standards
 * and going through the uniprot flat file database and extracting
 * all of the proteins on that list and then saving them in a
 * new file
 * @author Brian Risk
 *
 */

public class ExtractListFromDatabase {
	
	public static void main(String args[]) {
		U.p("Starting protein search...");
		File databaseFile = new File("tests/databases/uniprot_sprot.fasta");
		File listFile = new File("UPS/UPS-list.txt");
		File extractedProteinsFile = new File("UPS/extracted-proteins.txt");
		try {
			BufferedReader listReader = new BufferedReader(new FileReader(listFile));
			BufferedReader dbReader = new BufferedReader(new FileReader(databaseFile));
			PrintWriter extractedWriter = new PrintWriter(new BufferedWriter(new FileWriter(extractedProteinsFile)));
			
			//Extract the names of the proteins and store them in a list
			ArrayList<String> proteinNames = new ArrayList<String>();
			String proteinName = listReader.readLine();
			while (proteinName != null) {
				proteinName = proteinName.trim();
				//add on ">" symbol
				if (!proteinName.equals("")) proteinNames.add(">" + proteinName);
				proteinName = listReader.readLine();
			}
			
			//go through the database and find the members of this list and save them
			String dbLine = dbReader.readLine();
			boolean foundProtein = false;
			while (dbLine != null) {
				if (dbLine.startsWith(">")) {
					foundProtein = false;
					String name;
					for (int i = 0; i < proteinNames.size(); i++) {
						name = proteinNames.get(i);
						if (dbLine.startsWith(name)) {
							foundProtein = true;
							proteinNames.remove(i);
							extractedWriter.println(dbLine);
							break;
						}
					}
				} else {
					if (foundProtein) {
						extractedWriter.println(dbLine);
					}
				}
				dbLine = dbReader.readLine();
			}
			
			//List all unfound proteins
			if (proteinNames.size() > 0) {
				U.p("we did not find these proteins:");
				for (String name: proteinNames) {
					U.p(name);
				}
			}
			
			//close out our files
			listReader.close();
			dbReader.close();
			extractedWriter.flush();
			extractedWriter.close();
			
			U.p("Done.");
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
