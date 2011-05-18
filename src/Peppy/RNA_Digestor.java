package Peppy;

import java.util.ArrayList;
import java.util.Collections;

public class RNA_Digestor {
	
	RNA_Sequence rna;
	int length;
	int halfWindowSize = 60;
	int windowSize = 2 * halfWindowSize;
	int minimumIntronLength = 40;
	int maximumIntronLength = 50000;
	ArrayList<Peptide> peptides = new ArrayList<Peptide>();
	
	
	public RNA_Digestor(RNA_Sequence rna) {
		this.rna = rna;
		length = rna.getRNA_3to5().length;
		fullDigest();
	}
	
	
	public ArrayList<Peptide> getPeptides() {
		return peptides;
	}


	private ArrayList<Peptide> digest(
			int startIndex, 
			int stopIndex, 
			boolean isForward, 
			byte [] code
			) {
		
		return digest(startIndex, stopIndex, isForward, code, -1, rna.getStart(), length, null);
	}
	
	/**
	 * 
	 * @param startIndex
	 * @param stopIndex
	 * @param isForward
	 * @param codeChunk
	 * @param spliceLocation -1 if not a spliced sequence
	 * @param beginLocus
	 * @param intronLength
	 * @return
	 */
	private ArrayList<Peptide> digest(
			int startIndex, 
			int stopIndex, 
			boolean isForward, 
			byte [] codeChunk, 
			int spliceLocation, 
			int beginLocus, 
			int intronLength,
			int [] absoluteIndicies
			) {
		
		//init peptide list
		ArrayList<Peptide> peptides = new ArrayList<Peptide>();
		//a codon is a set of 3 nucleotide characters
		byte [] codon = new byte[3];
		//as we walk through the sequence, this index keeps track of which part of the codon array to fill
		int codonIndex = 0;
		//the amino acid is the translation of the codon
		char aminoAcid = 'x';
		char previousAminoAcid = 'x';
		//Where we store all of our forming peptides
		ArrayList<PeptideUnderConstruction> peptidesUnderConstruction = new ArrayList<PeptideUnderConstruction>();
		
		final int codonIncrement = 1;
		//NOTE commenting these as with RNA we have two separate forwards and reverse sequences
//		if (isForward) {codonIncrement = 1;}
//		else {codonIncrement = -1;}
		
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
				
				aminoAcid = Definitions.aminoAcidList[indexForCodonArray(codon, isForward)];
				
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

				//public Peptide(String acidSequence, int index, boolean forward, byte readingFrame, Sequence parentSequence, boolean isSpliced) {
				//Events which might start a new peptide
				if (aminoAcid == 'M') {
					inORF = true;
				} else {
					if (isBreak(aminoAcid)) {
						//add all peptides
						for (PeptideUnderConstruction puc: peptidesUnderConstruction) {
							//TODO add more info (start, stop, etc.) for peptide
							Peptide peptide = null;
							if (
									(spliceLocation == -1) ||
									(puc.getStartIndex() <= spliceLocation && acidIndex > spliceLocation)
									){
								//O, God!  keeping track of start and stop locations!  The madness consumes me!!!
								//all values initially set to -1 as that is not a valid index.  If you see a -1,
								//that means that value was never set -- which may be the case for the intron
								//indices.
								int peptideStartIndex = -1, peptideStopIndex = -1, intronStartIndex = -1, intronStopIndex = -1;
								if (spliceLocation == -1) {
									if (isForward) {
										peptideStartIndex = beginLocus + puc.getStartIndex();
										peptideStopIndex = beginLocus + acidIndex;
									} else {
										peptideStartIndex = beginLocus + codeChunk.length - puc.getStartIndex();
										peptideStopIndex = beginLocus + codeChunk.length - acidIndex;
									}	
								} else {
									peptideStartIndex = absoluteIndicies[puc.getStartIndex()];
									peptideStopIndex = absoluteIndicies[acidIndex];
									intronStartIndex = absoluteIndicies[spliceLocation];
									intronStopIndex = absoluteIndicies[spliceLocation + 1];
								}
								peptide = new Peptide(
										puc.getSequence(),
										peptideStartIndex,
										peptideStopIndex,
										intronStartIndex,
										intronStopIndex,
										isForward,
										null,
										(spliceLocation != -1)
										);
							}
							if (peptide != null) {
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
					}
					if (aminoAcid == '.') inORF = false;
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
				acidIndex += threeTimesCodonIncrement;
				previousAminoAcid = aminoAcid;
			}
			
			codon[codonIndex] = codeChunk[nucleotideIndex];	
			codonIndex++;
			nucleotideIndex += codonIncrement;	
		}	
		//adding all the remaining peptides under construction
		for (PeptideUnderConstruction puc: peptidesUnderConstruction) {
			Peptide peptide = null;
			if ((spliceLocation == -1) ||
					(puc.getStartIndex() <= spliceLocation && acidIndex > spliceLocation)){
				//TODO any way I can avoid exactly repeating the code from above?
				int peptideStartIndex = -1, peptideStopIndex = -1, intronStartIndex = -1, intronStopIndex = -1;
				if (spliceLocation == -1) {
					if (isForward) {
						peptideStartIndex = beginLocus + puc.getStartIndex();
						peptideStopIndex = beginLocus + acidIndex;
					} else {
						peptideStartIndex = beginLocus + codeChunk.length - puc.getStartIndex();
						peptideStopIndex = beginLocus + codeChunk.length - acidIndex;
					}	
				} else {
					peptideStartIndex = absoluteIndicies[puc.getStartIndex()];
					peptideStartIndex = absoluteIndicies[acidIndex];
				}
				peptide = new Peptide(
						puc.getSequence(),
						peptideStartIndex,
						peptideStopIndex,
						intronStartIndex,
						intronStopIndex,
						isForward,
						null,
						(spliceLocation != -1)
						);
			}
			if (peptide != null) {
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
		
		return peptides;
	}
	
	
	//does both normal and spliced digestion
	private void fullDigest() {
		//normal digestion
//		peptides.addAll(digest(0,length, true, rna.getRNA_5to3()));
//		peptides.addAll(digest(1,length, true, rna.getRNA_5to3()));
//		peptides.addAll(digest(2,length, true, rna.getRNA_5to3()));
//		peptides.addAll(digest(0,length, true, rna.getRNA_3to5()));
//		peptides.addAll(digest(1,length, true, rna.getRNA_3to5()));
//		peptides.addAll(digest(2,length, true, rna.getRNA_3to5()));
		
		//finding splices
		digestSplices(rna.getRNA_5to3(), true, rna.getForwardsStartLocations(), rna.getForwardsStopLocations());	
		digestSplices(rna.getRNA_3to5(), false, rna.getReverseStartLocations(), rna.getReverseStopLocations());
		
		Collections.sort(peptides);
		
	}
	
	private void digestSplices(byte [] code, boolean isForward, boolean[] starts, boolean[] stops) {
		byte [] splice = new byte[windowSize];
		int [] absoluteIndicies = new int[windowSize];
		for (int i = halfWindowSize; i < length - halfWindowSize; i++) {
			if (starts[i]) {
				for (int x = 0; x < halfWindowSize; x++) {
					splice[x] = code[i - halfWindowSize + x];
					if (isForward) {
						absoluteIndicies[x] = rna.getStart() + i - halfWindowSize + x;
					} else {
						absoluteIndicies[x] = rna.getStart() + length - (i - halfWindowSize + x);
					}
				}
				int jStop = i + maximumIntronLength;
				if (jStop > length - halfWindowSize) jStop = length - halfWindowSize;
				for (int j = i + minimumIntronLength; j < jStop; j++) {
					if (stops[j]) {
						for (int x = 0; x < halfWindowSize; x++) {
							splice[x + halfWindowSize] = code[j + x + 1];
							if (isForward) {
								absoluteIndicies[x + halfWindowSize] = rna.getStart() + j + x + 1;
							} else {
								absoluteIndicies[x + halfWindowSize] = rna.getStart() + length - (j + x + 1);
							}
						}
						//TODO starts and stops not handled correctly for 3to5
						peptides.addAll(digest(0, windowSize, isForward, splice, halfWindowSize - 1, rna.getStart() + i - halfWindowSize, j - i, absoluteIndicies));
						peptides.addAll(digest(1, windowSize, isForward, splice, halfWindowSize - 1, rna.getStart() + i - halfWindowSize, j - i, absoluteIndicies));
						peptides.addAll(digest(2, windowSize, isForward, splice, halfWindowSize - 1, rna.getStart() + i - halfWindowSize, j - i, absoluteIndicies));
					}
				}
			}
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
	private int indexForCodonArray(byte [] codon) {
		if (codon[0] == RNA_Sequence.BASE_A) {
			if (codon[1] == RNA_Sequence.BASE_A) {
				if (codon[2] == RNA_Sequence.BASE_A) {
					return 0;
				} else if (codon[2] == RNA_Sequence.BASE_C) {
					return 1;
				} else if (codon[2] == RNA_Sequence.BASE_G) {
					return 2;
				} else {
					return 3;
				}
			} else if (codon[1] == RNA_Sequence.BASE_C) {
				if (codon[2] == RNA_Sequence.BASE_A) {
					return 4;
				} else if (codon[2] == RNA_Sequence.BASE_C) {
					return 5;
				} else if (codon[2] == RNA_Sequence.BASE_G) {
					return 6;
				} else {
					return 7;
				}
			} else if (codon[1] == RNA_Sequence.BASE_G) {
				if (codon[2] == RNA_Sequence.BASE_A) {
					return 8;
				} else if (codon[2] == RNA_Sequence.BASE_C) {
					return 9;
				} else if (codon[2] == RNA_Sequence.BASE_G) {
					return 10;
				} else {
					return 11;
				}
			} else {
				if (codon[2] == RNA_Sequence.BASE_A) {
					return 12;
				} else if (codon[2] == RNA_Sequence.BASE_C) {
					return 13;
				} else if (codon[2] == RNA_Sequence.BASE_G) {
					return 14;
				} else {
					return 15;
				}
			}
		} else if (codon[0] == RNA_Sequence.BASE_C) {
			if (codon[1] == RNA_Sequence.BASE_A) {
				if (codon[2] == RNA_Sequence.BASE_A) {
					return 16;
				} else if (codon[2] == RNA_Sequence.BASE_C) {
					return 17;
				} else if (codon[2] == RNA_Sequence.BASE_G) {
					return 18;
				} else {
					return 19;
				}
			} else if (codon[1] == RNA_Sequence.BASE_C) {
				if (codon[2] == RNA_Sequence.BASE_A) {
					return 20;
				} else if (codon[2] == RNA_Sequence.BASE_C) {
					return 21;
				} else if (codon[2] == RNA_Sequence.BASE_G) {
					return 22;
				} else {
					return 23;
				}
			} else if (codon[1] == RNA_Sequence.BASE_G) {
				if (codon[2] == RNA_Sequence.BASE_A) {
					return 24;
				} else if (codon[2] == RNA_Sequence.BASE_C) {
					return 25;
				} else if (codon[2] == RNA_Sequence.BASE_G) {
					return 26;
				} else {
					return 27;
				}
			} else {
				if (codon[2] == RNA_Sequence.BASE_A) {
					return 28;
				} else if (codon[2] == RNA_Sequence.BASE_C) {
					return 29;
				} else if (codon[2] == RNA_Sequence.BASE_G) {
					return 30;
				} else {
					return 31;
				}
			}
		} else if (codon[0] == RNA_Sequence.BASE_G) {
			if (codon[1] == RNA_Sequence.BASE_A) {
				if (codon[2] == RNA_Sequence.BASE_A) {
					return 32;
				} else if (codon[2] == RNA_Sequence.BASE_C) {
					return 33;
				} else if (codon[2] == RNA_Sequence.BASE_G) {
					return 34;
				} else {
					return 35;
				}
			} else if (codon[1] == RNA_Sequence.BASE_C) {
				if (codon[2] == RNA_Sequence.BASE_A) {
					return 36;
				} else if (codon[2] == RNA_Sequence.BASE_C) {
					return 37;
				} else if (codon[2] == RNA_Sequence.BASE_G) {
					return 38;
				} else {
					return 39;
				}
			} else if (codon[1] == RNA_Sequence.BASE_G) {
				if (codon[2] == RNA_Sequence.BASE_A) {
					return 40;
				} else if (codon[2] == RNA_Sequence.BASE_C) {
					return 41;
				} else if (codon[2] == RNA_Sequence.BASE_G) {
					return 42;
				} else {
					return 43;
				}
			} else {
				if (codon[2] == RNA_Sequence.BASE_A) {
					return 44;
				} else if (codon[2] == RNA_Sequence.BASE_C) {
					return 45;
				} else if (codon[2] == RNA_Sequence.BASE_G) {
					return 46;
				} else {
					return 47;
				}
			}
		} else if (codon[0] == RNA_Sequence.BASE_U) {
			if (codon[1] == RNA_Sequence.BASE_A) {
				if (codon[2] == RNA_Sequence.BASE_A) {
					return 48;
				} else if (codon[2] == RNA_Sequence.BASE_C) {
					return 49;
				} else if (codon[2] == RNA_Sequence.BASE_G) {
					return 50;
				} else {
					return 51;
				}
			} else if (codon[1] == RNA_Sequence.BASE_C) {
				if (codon[2] == RNA_Sequence.BASE_A) {
					return 52;
				} else if (codon[2] == RNA_Sequence.BASE_C) {
					return 53;
				} else if (codon[2] == RNA_Sequence.BASE_G) {
					return 54;
				} else {
					return 55;
				}
			} else if (codon[1] == RNA_Sequence.BASE_G) {
				if (codon[2] == RNA_Sequence.BASE_A) {
					return 56;
				} else if (codon[2] == RNA_Sequence.BASE_C) {
					return 57;
				} else if (codon[2] == RNA_Sequence.BASE_G) {
					return 58;
				} else {
					return 59;
				}
			} else {
				if (codon[2] == RNA_Sequence.BASE_A) {
					return 60;
				} else if (codon[2] == RNA_Sequence.BASE_C) {
					return 61;
				} else if (codon[2] == RNA_Sequence.BASE_G) {
					return 62;
				} else {
					return 63;
				}
			}
		} else {
			return -1; //return STOP
		}
	}
	
	private int indexForCodonArray(byte [] codon, boolean forwards) {
		int out = indexForCodonArray(codon);
		if (out == -1) return 56; //if unknown, return STOP
		if (forwards) {
			return indexForCodonArray(codon);
		} else {
			return 63 - indexForCodonArray(codon);
		}
	}

	private boolean isBreak(char aminoAcid) {
		return (aminoAcid == '.' || aminoAcid == 'K' || aminoAcid == 'R');
	}

}
