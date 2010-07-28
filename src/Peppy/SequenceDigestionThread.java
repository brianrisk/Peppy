package Peppy;
import java.util.ArrayList;

import Utilities.U;

public class SequenceDigestionThread implements Runnable {
	
	ArrayList<Peptide> peptides = new ArrayList<Peptide>();
	NucleotideSequence nucleotideSequence;
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
	public SequenceDigestionThread(NucleotideSequence nucleotideSequence,
			byte frame, boolean forwards, int startIndex, int stopIndex) {
		this.nucleotideSequence = nucleotideSequence;
		this.frame = frame;
		this.forwards = forwards;
		this.startIndex = startIndex;
		this.stopIndex = stopIndex;
	}
	
	/**
	 * @param nucleotideSequence
	 * @param frame
	 * @param forwards
	 */
	public SequenceDigestionThread(NucleotideSequence nucleotideSequence,
			byte frame, boolean forwards) {
		this.nucleotideSequence = nucleotideSequence;
		this.frame = frame;
		this.forwards = forwards;
		this.startIndex = 0;
		this.stopIndex = nucleotideSequence.getSequence().length();
	}

	public ArrayList<Peptide> getPeptides( ) {
		if (peptides.size() == 0) {
			if (forwards) {
				digest(frame, nucleotideSequence.getSequence().length());
			} else {
				digest(nucleotideSequence.getSequence().length() - frame - 1, 0);	
			}
		}
		return peptides;
	}
	
	public void digest(int startIndex, int stopIndex) {
		//a codon is a set of 3 nucleotide characters
		char [] codon = new char[3];
		//as we walk through the sequence, this index keeps track of which part of the codon array to fill
		int codonIndex = 0;
		//the amino acid is the translation of the codon
		char aminoAcid = 'x';
		//Where we store all of our forming peptides
		ArrayList<PeptideUnderConstruction> peptidesUnderConstruction = new ArrayList<PeptideUnderConstruction>();
		
		final int codonIncrement;
		if (forwards) {codonIncrement = 1;}
		else {codonIncrement = -1;}
		
		int nucleotideIndex = startIndex;
		//acidIndex points to the beginning of where the acid is coded
		int acidIndex = startIndex;
		//so we don't have to keep doing this calculation
		int startIndexPlusThree = startIndex + (3 * codonIncrement);
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
				
				
				if (aminoAcid == 'M') {
					peptidesUnderConstruction.add(new PeptideUnderConstruction(acidIndex, 'M'));
					//sometimes the M is not added, so we're accounting for this here
					peptidesUnderConstruction.add(new PeptideUnderConstruction(acidIndex + (3 * codonIncrement)));
				} else {
					if (aminoAcid == '.' || aminoAcid == 'K' || aminoAcid == 'R') {
						peptidesUnderConstruction.add(new PeptideUnderConstruction(acidIndex + (3 * codonIncrement)));
						//add all peptides
						for (PeptideUnderConstruction puc: peptidesUnderConstruction) {
							Peptide peptide = new Peptide(puc.getSequence(), acidIndex, forwards, frame,  nucleotideSequence.getParentSequence());
							if (peptide.getMass() >= Properties.peptideMassThreshold) peptides.add(peptide);
						}
					}
				}
				
				//if stop, then clear out
				if (aminoAcid == '.') {
					peptidesUnderConstruction = new ArrayList<PeptideUnderConstruction>();
				}
				
				//remove all peptide under construction that have reached their maximum break count
				int size = peptidesUnderConstruction.size();
				for (int i = 0; i < size; i++) {
					PeptideUnderConstruction puc = peptidesUnderConstruction.get(i);
					if (puc.getBreakCount() > Properties.numberOfMissedCleavages) {
						peptidesUnderConstruction.remove(i);
						i--;
						size--;
					}
				}
				
				//our acid index can now move forward
				acidIndex += (3 * codonIncrement);
			}
			
			codon[codonIndex] = nucleotideSequence.getSequence().charAt(nucleotideIndex);	
			codonIndex++;
			nucleotideIndex += codonIncrement;	
		}	
		//adding all the remaining peptides under construction
		for (PeptideUnderConstruction puc: peptidesUnderConstruction) {
			Peptide peptide = new Peptide(puc.getSequence(), nucleotideIndex, forwards, frame,  nucleotideSequence.getParentSequence());
			if (peptide.getMass() >= Properties.peptideMassThreshold) peptides.add(peptide);
		}
	}
	
	/**
	 * The original digestion.  Can be deleted when you see fit.
	 * @param startIndex
	 * @param stopIndex
	 */
	public void digestORIGINAL(int startIndex, int stopIndex) {
		//a codon is a set of 3 nucleotide characters
		char [] codon = new char[3];
		//as we walk through the sequence, this index keeps track of which part of the codon array to fill
		int codonIndex = 0;
		//the amino acid is the translation of the codon
		char aminoAcid = 'x';
		//need to keep track of the last amino acid for trypsin reasons.
		//filling it with a nonsense value to start.
		char previousAminoAcid = 'x';
		//the index of our amino acid defined by the nucleotide codon
		int aminoAcidIndex;
		//peptide is the chain of amino acids we are building in the open reading frame
		StringBuffer proteinUnderConstruction = new StringBuffer();
		int peptideStartIndex = startIndex;
		String proteinString = "";
		

		
		//since we only want peptides in open reading frames (between a START and a STOP)
		boolean inOpenReadingFrame = false;
		
		final int codonIncrement;
		if (forwards) {codonIncrement = 1;}
		else {codonIncrement = -1;}
		
		
		//go through each character in the line, skipping ahead "frame" characters
		//this while loop is necessary (rather than a "for") to accommodate reading forwards/backwards
		int nucleotideIndex = startIndex;
		while (nucleotideIndex != stopIndex) {
			
			
			//if codonIndex is 3 that means our codon is full and we can determine the amino acid
			if (codonIndex == 3) {
				codonIndex = 0;
				aminoAcidIndex = indexForCodonArray(codon, forwards);
				aminoAcid = Definitions.aminoAcidList[aminoAcidIndex];
				
				//determine if we're in an open reading frame
				
				//if previous amino is stop, we may want to dump the created protein into the protein digester
				if (previousAminoAcid == '.') {
					if (inOpenReadingFrame || !Properties.onlyUsePeptidesInOpenReadingFrames) {
						proteinString = proteinUnderConstruction.toString();
						peptides.addAll(ProteinDigestion.getPeptidesFromProteinString(proteinString, peptideStartIndex, forwards, frame, nucleotideSequence.getParentSequence()));
					}
					inOpenReadingFrame = false;
				}
				if (previousAminoAcid == 'M' && !inOpenReadingFrame) {
					inOpenReadingFrame = true;
					proteinUnderConstruction = new StringBuffer();
					peptideStartIndex = nucleotideIndex;
					//shifting to correspond with GFS
					if (!forwards) {peptideStartIndex += 1;}
				}
				
				if (inOpenReadingFrame || !Properties.onlyUsePeptidesInOpenReadingFrames) {proteinUnderConstruction.append(aminoAcid);}
				//peptideMass += Definitions.aminoAcidMassesAverage[aminoAcidIndex];
				previousAminoAcid = aminoAcid;
			}
			
//			proteinString = proteinUnderConstruction.toString();
//			peptides.addAll(ProteinDigestion.getPeptidesFromProteinString(proteinString, peptideStartIndex, forwards, frame, nucleotideSequence.getParentSequence()));
			
			codon[codonIndex] = nucleotideSequence.getSequence().charAt(nucleotideIndex);	
			codonIndex++;
			nucleotideIndex += codonIncrement;	
		}
		
		
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


