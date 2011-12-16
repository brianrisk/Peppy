package Peppy;
import java.util.ArrayList;

public class DigestionThread_DNA implements Runnable {
	
	ArrayList<Peptide> peptides = new ArrayList<Peptide>();
	Nucleotides nucleotideSequence;
	int frame;
	boolean isForwardsStrand;
	int startIndex;
	int stopIndex;
	boolean reverseDatabase;

	/**
	 * @param nucleotideSequence
	 * @param frame
	 * @param forwardsStrand
	 */
	public DigestionThread_DNA(Nucleotides nucleotideSequence,
			int frame, boolean forwardsStrand, int startIndex, int stopIndex, boolean reverseDatabase) {
		
		/* because of the digestion window overlap, sometimes this will be < 0 */
		if (startIndex < 0) startIndex = 0;
		
		/* set our variables */
		this.nucleotideSequence = nucleotideSequence;
		this.frame = frame;
		this.isForwardsStrand = forwardsStrand;
		this.startIndex = startIndex;
		this.stopIndex = stopIndex;
		this.reverseDatabase = reverseDatabase;
	}


	public void run() {
		digest();
	}
	

	public ArrayList<Peptide> digest() {
		/* go through each character in the line, skipping ahead "frame" characters */
		if (isForwardsStrand) {
			 return digest(startIndex + frame, stopIndex);
		} else {
			/* the "startIndex -1" is because we go while our 
			 * index does not equal to that number.  Therefore it goes
			 * right up to startIndex  */
			return digest(stopIndex - frame - 1, startIndex - 1);	
		}
	}
	
	
	private ArrayList<Peptide> digest(int startPosition, int stopPosition) {
		ArrayList<Protein> proteins = translateToProteins(startPosition, stopPosition);
		peptides = Sequence_Protein.getPeptidesFromListOfProteins(proteins);
		return peptides;
	}
	

	/**
	 * creates an ArrayList of Proteins where STOPs mark the end
	 * of each protein
	 * @param startPosition
	 * @param stopPosition
	 */
	private ArrayList<Protein> translateToProteins(int startPosition, int stopPosition) {
		ArrayList<Protein> proteins = new ArrayList<Protein>();
		char [] codon = new char[3];
		char aminoAcid;
		int mod = 0;
		StringBuffer buildingProtein = new StringBuffer();
		Sequence_DNA sequence_DNA = nucleotideSequence.getParentSequence();
		String name = nucleotideSequence.getSequenceDescription();
		int increment = 1;
		if (!isForwardsStrand) increment = -1;
		int index;
		int proteinStart = startPosition;
		for (index = startPosition; index != stopPosition; index += increment) {
			codon[mod] = nucleotideSequence.getSequence().charAt(index);
			if (mod == 2) {
				aminoAcid = Definitions.aminoAcidList[indexForCodonArray(codon, isForwardsStrand)];
				buildingProtein.append(aminoAcid);
				if (aminoAcid == '.') {
					if (buildingProtein.length() > 3) {
						if (reverseDatabase) buildingProtein.reverse();
						proteins.add(new Protein(name, proteinStart, buildingProtein.toString(), false, -1, -1, isForwardsStrand, sequence_DNA));
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
			proteins.add(new Protein(name, proteinStart, buildingProtein.toString(), false, -1, -1, isForwardsStrand, sequence_DNA));
		}
		return proteins;
	}


	public static int indexForCodonArray(char [] codon, boolean forwards) {
		int out = indexForCodonArray(codon);
		if (out == -1) return 56; //if unknown, return STOP
		if (forwards) {
			return indexForCodonArray(codon);
		} else {
			return 63 - indexForCodonArray(codon);
		}
	}
	
	/*
	 * TODO this method needs to make sure that if any unknown characters are found then 
	 * STOP (56) is returned.  Right now that only happens with the last character.
	 */
	/**
	 * Same as other, but assumes that the direction is "forwards"
	 * @param codon
	 * @return
	 */
	public static int indexForCodonArray(char [] codon) {
		if (codon[0] == 'A') {
			if (codon[1] == 'A') {
				if (codon[2] == 'A') {
					return 0;
				} else if (codon[2] == 'C') {
					return 1;
				} else if (codon[2] == 'G') {
					return 2;
				} else {
					return 3;
				}
			} else if (codon[1] == 'C') {
				if (codon[2] == 'A') {
					return 4;
				} else if (codon[2] == 'C') {
					return 5;
				} else if (codon[2] == 'G') {
					return 6;
				} else {
					return 7;
				}
			} else if (codon[1] == 'G') {
				if (codon[2] == 'A') {
					return 8;
				} else if (codon[2] == 'C') {
					return 9;
				} else if (codon[2] == 'G') {
					return 10;
				} else {
					return 11;
				}
			} else {
				if (codon[2] == 'A') {
					return 12;
				} else if (codon[2] == 'C') {
					return 13;
				} else if (codon[2] == 'G') {
					return 14;
				} else {
					return 15;
				}
			}
		} else if (codon[0] == 'C') {
			if (codon[1] == 'A') {
				if (codon[2] == 'A') {
					return 16;
				} else if (codon[2] == 'C') {
					return 17;
				} else if (codon[2] == 'G') {
					return 18;
				} else {
					return 19;
				}
			} else if (codon[1] == 'C') {
				if (codon[2] == 'A') {
					return 20;
				} else if (codon[2] == 'C') {
					return 21;
				} else if (codon[2] == 'G') {
					return 22;
				} else {
					return 23;
				}
			} else if (codon[1] == 'G') {
				if (codon[2] == 'A') {
					return 24;
				} else if (codon[2] == 'C') {
					return 25;
				} else if (codon[2] == 'G') {
					return 26;
				} else {
					return 27;
				}
			} else {
				if (codon[2] == 'A') {
					return 28;
				} else if (codon[2] == 'C') {
					return 29;
				} else if (codon[2] == 'G') {
					return 30;
				} else {
					return 31;
				}
			}
		} else if (codon[0] == 'G') {
			if (codon[1] == 'A') {
				if (codon[2] == 'A') {
					return 32;
				} else if (codon[2] == 'C') {
					return 33;
				} else if (codon[2] == 'G') {
					return 34;
				} else {
					return 35;
				}
			} else if (codon[1] == 'C') {
				if (codon[2] == 'A') {
					return 36;
				} else if (codon[2] == 'C') {
					return 37;
				} else if (codon[2] == 'G') {
					return 38;
				} else {
					return 39;
				}
			} else if (codon[1] == 'G') {
				if (codon[2] == 'A') {
					return 40;
				} else if (codon[2] == 'C') {
					return 41;
				} else if (codon[2] == 'G') {
					return 42;
				} else {
					return 43;
				}
			} else {
				if (codon[2] == 'A') {
					return 44;
				} else if (codon[2] == 'C') {
					return 45;
				} else if (codon[2] == 'G') {
					return 46;
				} else {
					return 47;
				}
			}
		} else if (codon[0] == 'T') {
			if (codon[1] == 'A') {
				if (codon[2] == 'A') {
					return 48;
				} else if (codon[2] == 'C') {
					return 49;
				} else if (codon[2] == 'G') {
					return 50;
				} else {
					return 51;
				}
			} else if (codon[1] == 'C') {
				if (codon[2] == 'A') {
					return 52;
				} else if (codon[2] == 'C') {
					return 53;
				} else if (codon[2] == 'G') {
					return 54;
				} else {
					return 55;
				}
			} else if (codon[1] == 'G') {
				if (codon[2] == 'A') {
					return 56;
				} else if (codon[2] == 'C') {
					return 57;
				} else if (codon[2] == 'G') {
					return 58;
				} else {
					return 59;
				}
			} else {
				if (codon[2] == 'A') {
					return 60;
				} else if (codon[2] == 'C') {
					return 61;
				} else if (codon[2] == 'G') {
					return 62;
				} else {
					return 63;
				}
			}
		} else {
			return -1; //return STOP
		}
	}


	public ArrayList<Peptide> getPeptides( ) {
		return peptides;
	}

}


