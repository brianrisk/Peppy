package Peppy;

import java.util.ArrayList;

public class ShortNucleotideDigestionThread  implements Runnable {
	
	private static final boolean [] strandValues = {true,false};
	private static final int [] frameValues = {0,1,2};
	
	NucleotideSequence nucleotideSequence;
	ShortNucleotideDigestionServer shortNucleotideDigestionServer;
	boolean reverseDatabase;
	
	/**
	 * @param proteins
	 * @param spectrum
	 */
	public ShortNucleotideDigestionThread(NucleotideSequence nucleotideSequence, ShortNucleotideDigestionServer shortNucleotideDigestionServer, boolean reverseDatabase) {
		this.nucleotideSequence = nucleotideSequence;
		this.shortNucleotideDigestionServer = shortNucleotideDigestionServer;
		this.reverseDatabase = reverseDatabase;
	}
	

	public void run() {
		while (nucleotideSequence != null) {
			//return results, get new task
			shortNucleotideDigestionServer.report(translate());
			nucleotideSequence = shortNucleotideDigestionServer.getNextElementToProcess();
		}
	}

	public ArrayList<Protein> translate() {
		ArrayList<Protein> proteins = new ArrayList<Protein>();
		
		/* getting out if the sequence isn't long enough for even one AA */
		if (nucleotideSequence.getSequence().length() <3) return proteins;
		
		char [] codon = new char[3];
		char aminoAcid;
		int mod;
		StringBuffer buildingProtein;
		Sequence_DNA sequence_DNA = nucleotideSequence.getParentSequence();
		String name = nucleotideSequence.getSequenceDescription();
		int index, increment, startIndex, stopIndex;
		for (boolean isForwardsStrand: strandValues) {
			for (int frameValue: frameValues) {
				mod = 0;
				if (isForwardsStrand) {
					increment = 1;
					startIndex = frameValue;
					stopIndex = nucleotideSequence.getSequence().length();
				} else {
					increment = -1;
					startIndex = nucleotideSequence.getSequence().length() - frameValue - 1;
					stopIndex = -1;
				}
				int proteinStart = 0;
				
				/* init the protein we are building */
				buildingProtein = new StringBuffer();
				
				for (index = startIndex; index != stopIndex; index += increment) {
					codon[mod] = nucleotideSequence.getSequence().charAt(index);
					if (mod == 2) {
						aminoAcid = Definitions.aminoAcidList[DigestionThread_DNA.indexForCodonArray(codon, isForwardsStrand)];
						buildingProtein.append(aminoAcid);
						if (aminoAcid == '.') {
							if (buildingProtein.length() > 3) {
								if (reverseDatabase) {
									/* remove the '.' that was just added above */
									buildingProtein.deleteCharAt(buildingProtein.length() - 1);
									
									/* reverse the string! */
									buildingProtein.reverse();
								}
								proteins.add(new Protein(name, proteinStart, buildingProtein.toString(), false, -1, -1, isForwardsStrand, sequence_DNA, reverseDatabase));
							}
							buildingProtein = new StringBuffer();
							proteinStart = index + increment;
						}
						
						/* reset mod */
						mod = 0;
					} else {
						mod++;
					}
				}
				
				if (buildingProtein.length() > 3) {
					proteins.add(new Protein(name, proteinStart, buildingProtein.toString(), false, -1, -1, isForwardsStrand, sequence_DNA, reverseDatabase));
				}
			}
		}
		return proteins;
	}

	
	

}
