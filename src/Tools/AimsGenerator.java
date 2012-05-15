package Tools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

import Peppy.Peptide;
import Peppy.Sequence;
import Peppy.Sequence_Protein;
import Peppy.U;

public class AimsGenerator {
	
	public static void main(String[] args) {
		
		/* get the target sequence  */
		U.p("Give the full file path to the target protein:");
		String targetFilePath = U.in();
		File targetFile = new File(targetFilePath);
		if (!targetFile.exists()) {
			U.p("invalid file path.  Exiting.");
			System.exit(1);
		}
		Sequence tartetSequence = new Sequence_Protein(targetFile);
		
		/* get the reference sequence */
		U.p("\rGive the full file path to the reference protein");
		String referenceFilePath = U.in();
		File referenceFile = new File(referenceFilePath);
		if (!referenceFile.exists()) {
			U.p("invalid file path.  Exiting.");
			System.exit(1);
		}
		Sequence referenceSequence = new Sequence_Protein(referenceFile);
		
		/* get our target peptides */
		U.p("\rloading proteins");
		Hashtable<String, Peptide> targetPeptides = getPeptideStringsFromSequence(tartetSequence);
		ArrayList<String> targetStrings = new ArrayList<String>(targetPeptides.keySet());
		Hashtable<String, Peptide> referencePeptides = getPeptideStringsFromSequence(referenceSequence);
		
		
		/* collect our unique peptides for AIMS */
		ArrayList<Peptide> aimsPeptides = new ArrayList<Peptide>(targetPeptides.size());
		for (String peptideSequence: targetStrings) {
			if (referencePeptides.get(peptideSequence) == null) {
				aimsPeptides.add(targetPeptides.get(peptideSequence));
			}
		}
		
		
		/* save our report */
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(targetFile.getName() +  " AIMS.txt")));
			Collections.sort(aimsPeptides);
			for (Peptide peptide: aimsPeptides) {
				String line = peptide.getAcidSequenceString();
				line += "\t" + peptide.getMass();
				line += "\t" + peptide.getProtein().getName();
				pw.println(line);
			}
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		U.p("\rdone");
	}
	
	private static Hashtable<String, Peptide> getPeptideStringsFromSequence(Sequence sequence) {
		ArrayList<Peptide> targetPeptides =  sequence.extractAllPeptides(false);
		Hashtable<String, Peptide> targetHash = new Hashtable<String, Peptide>();
		for (Peptide peptide: targetPeptides) {
			targetHash.put(peptide.getAcidSequenceString(), peptide);
		}
		return targetHash;
	}
	

}
