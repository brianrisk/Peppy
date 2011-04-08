package Peppy;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;

import Utilities.U;

/**
 * Meant to be a tool to output FASTA and CSV 
 * databases for predicted alternatives splicings 
 * on a chromosome region.
 * @author Brian Risk
 *
 */
public class AltSpliceTool {
	
	static Sequence sequence;
	static int startIndex;
	static int stopIndex;
	
	/**
	 * args[] contains three parameters:
	 * 1) The path to the chromosome
	 * 2) the start index
	 * 3) the stop index
	 * 
	 * This method will then find all probable alternative
	 *  splicings along that region and output both a CSV 
	 *  and FASTA file
	 * @param args
	 */
	public static void main(String args[]) {
		//make sure we have the correct number of arguments
		if (args.length != 3) {
			U.p("There must be three arguments:  The file path to the chromosome; the start index; the stop index");
			System.exit(1);
		}
		
//		if (args.length == 0) {
//			String args2[] = {"/Users/risk2/PeppyOverflow/chr12.fa",  "65160417",  "65210417"};
//			args = args2;
//		}
		
		//make sure this file exists
		File chromosomeFile = new File(args[0]);
		if (!chromosomeFile.exists()) {
			U.p("No file exists at location: " + args[0]);
			System.exit(1);
		}
		
		//make sure start location is an integer
		Integer startLocationInteger = null;
		try {
			startLocationInteger = new Integer(args[1]);
		} catch (NumberFormatException nfe) {
			U.p("The value entered for the start location, \"" + args[1] + "\" is not an ingeter."); 
		}
		
		//make sure stop location is an integer
		Integer stopLocationInteger = null;
		try {
			stopLocationInteger = new Integer(args[2]);
		} catch (NumberFormatException nfe) {
			U.p("The value entered for the stop location, \"" + args[2] + "\" is not an ingeter."); 
		}
		
		//Now that the values have been checked, s up our variables
		ArrayList<Sequence> sequences = Sequence.loadSequences(chromosomeFile);
		sequence = sequences.get(0);
		startIndex = startLocationInteger.intValue();
		stopIndex = stopLocationInteger.intValue();
		String chrName = U.getFileNameWithoutSuffix(chromosomeFile);
		
		//performing alt splicing
		U.p("digesting " + chrName + " from " + startIndex + " to " + stopIndex + "...");
		RNA_Sequence rna = new RNA_Sequence(sequence.getNucleotideSequences().get(0), startIndex, stopIndex);
		rna.printStats();
		RNA_Digestor rnaDigestor = new RNA_Digestor(rna);
		ArrayList<Peptide> peptides  = rnaDigestor.getPeptides();
		U.p("peptide tally: " + peptides.size());
		
		//save the peptides to FASTA
		U.p("saving FASTA file...");
		File peptideFASTAFile = new File(chrName + "_" + startIndex + "-" + stopIndex + ".fasta");
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(peptideFASTAFile)));
			for (Peptide peptide: peptides) {
				if (peptide.isSpliced()) {
					pw.println(
							">" + chromosomeFile.getName() + 
							"," + peptide.isForward() + 
							"," + peptide.isSpliced() + 
							"," + peptide.getStartIndex() + 
							"," + peptide.getStopIndex() + 
							"," + peptide.getIntronStartIndex() + 
							"," + peptide.getIntronStopIndex()
							);
					pw.println(peptide.getAcidSequenceString());
				}
			}
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//save the peptides and masses to CSV
		U.p("saving CSV file...");
		File peptideCSVFile = new File(chrName + "_" + startIndex + "-" + stopIndex + ".csv");
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(peptideCSVFile)));
			double charge1, charge2, charge3, charge4, charge5;
			for (Peptide peptide: peptides) {
				if (peptide.isSpliced()) {
					charge1 = peptide.getMass() + Definitions.HYDROGEN_MONO;
					charge2 = peptide.getMass() / 2 + Definitions.HYDROGEN_MONO;
					charge3 = peptide.getMass() / 3 + Definitions.HYDROGEN_MONO;
					charge4 = peptide.getMass() / 4 + Definitions.HYDROGEN_MONO;
					charge5 = peptide.getMass() / 5 + Definitions.HYDROGEN_MONO;
					pw.println(
							peptide.getAcidSequenceString() + "\t" + 
							charge1 + "\t" + 
							charge2 + "\t" + 
							charge3 + "\t" + 
							charge4 + "\t" + 
							charge5
							);
				}
			}
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//done!
		U.p("done!");
		
	}

}
