package Peppy;
import java.util.ArrayList;

public class DNA_DigestionThread implements Runnable {
	
	ArrayList<Peptide> peptides = new ArrayList<Peptide>();
	DNA_Sequence nucleotideSequence;
	byte frame;
	boolean forwards;
	int startIndex;
	int stopIndex;

	public void run() {
		//go through each character in the line, skipping ahead "frame" characters
		if (forwards) {
			digest(startIndex + frame, stopIndex);
		} else {
			digest(stopIndex - frame - 1, startIndex);	
		}
	}
	

	/**
	 * @param nucleotideSequence
	 * @param frame
	 * @param forwards
	 */
	public DNA_DigestionThread(DNA_Sequence nucleotideSequence,
			byte frame, boolean forwards, int startIndex, int stopIndex) {
		if (startIndex < 0) startIndex = 0;
		this.nucleotideSequence = nucleotideSequence;
		this.frame = frame;
		this.forwards = forwards;
		this.startIndex = startIndex;
		this.stopIndex = stopIndex;
	}
	
	public ArrayList<Peptide> getPeptides( ) {
		return peptides;
	}
	
	
	public void digest(int startIndex, int stopIndex) {
		//a codon is a set of 3 nucleotide characters
		char [] codon = new char[3];
		//as we walk through the sequence, this index keeps track of which part of the codon array to fill
		int codonIndex = 0;
		//the amino acid is the translation of the codon
		char aminoAcid = 'x';
		char previousAminoAcid = 'x';
		//Where we store all of our forming peptides
		ArrayList<PeptideUnderConstruction> peptidesUnderConstruction = new ArrayList<PeptideUnderConstruction>();
		
		final int codonIncrement;
		if (forwards) {codonIncrement = 1;}
		else {codonIncrement = -1;}
		
		int threeTimesCodonIncrement = (3 * codonIncrement);
		
		//account for digestionFrameOverlap making startIndex < 0
		if (startIndex < 0) startIndex = 0;
		int nucleotideIndex = startIndex;
		//acidIndex points to the beginning of where the acid is coded
		int acidIndex = startIndex;
		//so we don't have to keep doing this calculation
		int startIndexPlusThree = startIndex + threeTimesCodonIncrement;
		
		//if we are interested in ORFs (Open Reading Frames)
		boolean inORF = false;
		
		
		while (nucleotideIndex != stopIndex) {
			
			
			//if codonIndex is 3 that means our codon is full and we can determine the amino acid
			if (codonIndex == 3) {
				codonIndex = 0;
				
				aminoAcid = Definitions.aminoAcidList[indexForCodonArray(codon, forwards)];
				
				//if this is the beginning, start a peptide
				if (nucleotideIndex == startIndexPlusThree) {
					peptidesUnderConstruction.add(new PeptideUnderConstruction(acidIndex, aminoAcid));
				} else {
					//add this amino acid to all peptides under construction
					for (PeptideUnderConstruction puc: peptidesUnderConstruction) {
						puc.addAminoAcid(aminoAcid);
					}
				}
				
				if ((aminoAcid == 'M') ||  // start a new peptide at M
					(previousAminoAcid == 'M' && aminoAcid != 'M') || // handle possible N-terminal methionine truncation products
					(isBreak(previousAminoAcid) && aminoAcid != 'M'))  // Create new peptides after a break, but only if we wouldn't have created a new one with M already
				{		
					peptidesUnderConstruction.add(new PeptideUnderConstruction(acidIndex, aminoAcid));
				}

				
				//Events which might start a new peptide
				if (aminoAcid == 'M') {
					inORF = true;
				} else {
					if (isBreak(aminoAcid)) {
						//if the next codon is not for STOP then add all peptides
						for (PeptideUnderConstruction puc: peptidesUnderConstruction) {
							Peptide peptide = new Peptide(puc.getSequence(), puc.getStartIndex(), forwards,  nucleotideSequence.getParentSequence());
							if (peptide.getMass() >= Properties.peptideMassThreshold) {
								if (Properties.onlyUsePeptidesInOpenReadingFrames) {
									if (inORF) {
										peptides.add(peptide);
									}
								} else {
									peptides.add(peptide);
								}
							}
						}
					}
					if (aminoAcid == '.') inORF = false;
				}
				
				//if stop, then clear out
				if (aminoAcid == '.') {
					peptidesUnderConstruction = new ArrayList<PeptideUnderConstruction>();
				}
				
				//remove all peptide under construction that have reached their maximum break count
				int size = peptidesUnderConstruction.size();
				for (int pucIndex = 0; pucIndex < size; pucIndex++) {
					PeptideUnderConstruction puc = peptidesUnderConstruction.get(pucIndex);
					if (puc.getBreakCount() > Properties.numberOfMissedCleavages) {
						peptidesUnderConstruction.remove(pucIndex);
						pucIndex--;
						size--;
					}
				}
				
				//our acid index can now move forward
				acidIndex += threeTimesCodonIncrement;
				previousAminoAcid = aminoAcid;
			}
			
			codon[codonIndex] = nucleotideSequence.getSequence().charAt(nucleotideIndex);	
			codonIndex++;
			nucleotideIndex += codonIncrement;	
		}	
		//adding all the remaining peptides under construction
		for (PeptideUnderConstruction puc: peptidesUnderConstruction) {
			Peptide peptide = new Peptide(puc.getSequence(), puc.getStartIndex(), forwards,  nucleotideSequence.getParentSequence());
			if (peptide.getMass() >= Properties.peptideMassThreshold) {
				if (Properties.onlyUsePeptidesInOpenReadingFrames) {
					if (inORF) {
						peptides.add(peptide);
					}
				} else {
					peptides.add(peptide);
				}
			}
		}
	}
	
	private boolean isBreak(char aminoAcid) {
		return (aminoAcid == '.' || aminoAcid == 'K' || aminoAcid == 'R');
	}
	
	
	


	private int indexForCodonArray(char [] codon, boolean forwards) {
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
	private int indexForCodonArray(char [] codon) {
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

}


