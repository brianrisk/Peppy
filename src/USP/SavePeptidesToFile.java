package USP;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;

import Peppy.Peptide;
import Peppy.ProteinDigestion;
import Utilities.U;

/**
 * takes our file of extracted proteins and then
 * saves a tab delimited list of sorted peptides and their properties
 * @author Brian Risk
 *
 */
public class SavePeptidesToFile {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		U.p("Saving peptides to file");
		ArrayList<Peptide> peptides = ProteinDigestion.getPeptidesFromProteinFile(new File("USP/extracted-proteins.txt"));
		Collections.sort(peptides);
		try {
			//printing full peptide list
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("USP/peptide-list.txt")));
			for (Peptide peptide: peptides) {
				pw.println(peptide.getAcidSequence());
			}
			
			//closing our file
			pw.flush();
			pw.close();
			
			//remove ununique massed peptides
			int listSize = peptides.size();
			double previousMass = 0;
			boolean justChanged = true;
			PrintWriter redundantMassList = new PrintWriter(new BufferedWriter(new FileWriter("USP/redundant-mass-list.txt")));
			for (int i = 0; i < listSize; i++) {
				Peptide pep = peptides.get(i);
				if (pep.getMass() - previousMass < 0.01) {
					if (justChanged) {
						redundantMassList.println(pep.getMass());
						justChanged = false;
						peptides.remove(i-1);
						i--;
						listSize--;
					}
					peptides.remove(i);
					i--;
					listSize--;
				} else {
					previousMass = pep.getMass();
					justChanged = true;
				}
			}
			
			//printing unique mass peptide list
			pw = new PrintWriter(new BufferedWriter(new FileWriter("USP/peptide-unique-mass-list.txt")));
			for (Peptide peptide: peptides) {
				pw.println(peptide.getAcidSequence());
			}
			
			//closing our file
			redundantMassList.flush();
			redundantMassList.close();
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		U.p("done");
	}

}
