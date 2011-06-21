package Tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import Utilities.U;

public class ConvertDATtoFASTA {
	
	public static void main(String args[]) {
		convertDATtoFASTA(new File(args[0]));
	}
	
	/**
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
