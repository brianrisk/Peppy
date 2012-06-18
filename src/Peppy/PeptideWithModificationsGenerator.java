package Peppy;

import java.util.ArrayList;

/**
 * For any given peptide there may exist a host of potentially modified form.
 * 
 * This class takes a peptide and a set of variable modifications
 * 
 * Copyright 2012, Brian Risk
 * Released under the Netscape Public License
 * 
 * 
 * @author Brian Risk
 *
 */
public class PeptideWithModificationsGenerator {
	
	Peptide peptide;
	ArrayList<ModificationVariable> modificaitons;
	
	
	/**
	 * see if it works
	 * @param args
	 */
	public static void main(String args[]) {
		Peptide peptide = new Peptide("MGGMSSK");
		ArrayList<ModificationVariable> mods = new ArrayList<ModificationVariable>();
		mods.add(new ModificationVariable(AminoAcids.M, 16));
		mods.add(new ModificationVariable(AminoAcids.M, 7));
		mods.add(new ModificationVariable(AminoAcids.K, 5));
		
		PeptideWithModificationsGenerator pmg = new PeptideWithModificationsGenerator(peptide, mods);
		ArrayList<PeptideWithModifications> peptides = pmg.getModificaitonVariaitons();
//		Collections.sort(peptides);
		for (PeptideWithModifications modPep: peptides) {
			double [] modArray = modPep.getModifications();
			StringBuffer toPrint = new StringBuffer();
			for (double val: modArray) {
				toPrint.append(val);
				toPrint.append(',');
			}
			U.p(toPrint);
		}
	}
	
	
	public PeptideWithModificationsGenerator(Peptide peptide,
			ArrayList<ModificationVariable> modificaitons) {
		this.peptide = peptide;
		this.modificaitons = modificaitons;
	}
	
	public ArrayList<PeptideWithModifications> getModificaitonVariaitons() {
		ArrayList<PeptideWithModifications> out = new ArrayList<PeptideWithModifications> ();
		byte [] sequence = peptide.getAcidSequence();
		String peptideString = peptide.getAcidSequenceString();
		
		int [] modificationIndices = new int [sequence.length];
		for (int index = 0; index < modificationIndices.length; index++) {
			modificationIndices[index] = -1;
		}
		int residueIndex = getNextRadixArray(sequence, 0, modificationIndices);
		
		while (residueIndex != -1) {
			double [] modificaitonArray = new double [sequence.length];
			for (int index = 0; index < modificaitonArray.length; index++) {
				if (modificationIndices[index] == -1 ) {
					modificaitonArray[index] = 0;
				} else {
					if (sequence[index] == modificaitons.get(modificationIndices[index]).getAminoAcid()) {
						modificaitonArray[index] = modificaitons.get(modificationIndices[index]).getMass();
					}
				}
			}
			out.add(new PeptideWithModifications(peptideString, modificaitonArray));
			

			residueIndex = getNextRadixArray(sequence, 0, modificationIndices);
		}
		

		return out;
	}
	
	
	
	private int getNextRadixArray(byte [] sequence, int residueIndex, int [] modificationIndices) {
		
		/* keep going until we find mod configuration that actually applies to the amino acids*/
		boolean keepGoing = true;
		
		/* if we went from one digit to the next in this step */
		boolean carryTheZero = false;
		
		while (keepGoing) {
			/* increment the indicies being sure to carry the zero */
			modificationIndices[residueIndex]++;
			
			/* e.g. 9999 to 10000 loop would repeat 4 times */
			while (modificationIndices[residueIndex] == modificaitons.size()) {
				modificationIndices[residueIndex] = -1;
				residueIndex++;
				carryTheZero = true;
				
				/* exit if we have reached the very end */
				if (residueIndex == modificationIndices.length) {
					return -1;
				}
			}
			
			/* if we went from one digit to the next, we increment the most sig.  e.g. 1399 to 1400 */
			if (carryTheZero) {
				modificationIndices[residueIndex]++;
			}
			
			/* we have reached the end */
			if (modificationIndices[residueIndex] == modificaitons.size()) {
				return -1;
			}
			
			/* see if we should keep going */
			if (modificationIndices[residueIndex] == -1) {
				keepGoing = true;
			} else {
				/* keep going if this mod has nothing to do with this amino acid */
				keepGoing = (sequence[residueIndex] != modificaitons.get(modificationIndices[residueIndex]).getAminoAcid());
			}
			

		}
		
		/* if we went from one digit to the next, we've got to go back to the least significant.  999 goes to 1000, but then that's back to 1001,1002 and so on */
		if (carryTheZero) {
			residueIndex = 0;
		}
		
		
		
		return residueIndex;
	}
	
	

}
