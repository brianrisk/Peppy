package Tools;

import Peppy.*;

import java.io.*;
import java.util.ArrayList;

/**
 * Save a list of digested peptides
 * @author Brian Risk
 *
 */
public class ExportPeptideList {
	
	public static void main(String args[]) {
		init("properties.txt");
		
		U.p("exporting peptide list");
		//Get references to our sequence files -- no nucleotide data is loaded at this point
		ArrayList<Sequence> sequences = Sequence.loadSequenceFiles(Properties.sequenceDirectoryOrFile);
		
		try {
			//loop through each sequence in the sequences ArrayList
			File peptidesFolder = new File ("exported peptides");
			peptidesFolder.mkdir();
			for (Sequence sequence: sequences) {
				U.p("working on sequence " + sequence.getSequenceFile().getName());
				PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(peptidesFolder, sequence.getSequenceFile().getName()))));
				
				ArrayList<Peptide> peptides = null;
				peptides = sequence.extractMorePeptides(false);
				while (peptides != null) {
					for (Peptide peptide: peptides) {
						pw.println(peptide);
					}
					peptides = sequence.extractMorePeptides(false);
				}
				peptides = null;
				System.gc();
				pw.flush();
				pw.close();
				U.p(sequence.getSequenceFile().getName() + " digested.");
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	public static void init(String propertiesFile) {
		System.setProperty("java.awt.headless", "true"); 
		Properties.loadProperties(propertiesFile);
		AminoAcids.init();
	}

}
