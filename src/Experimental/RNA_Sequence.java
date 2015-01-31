package Experimental;

import Peppy.NucleotideSequence;
import Peppy.SequenceNucleotide;
import Peppy.U;

/**
 * 
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class RNA_Sequence {
	
	SequenceNucleotide sequenceFile;
	NucleotideSequence DNA;
	byte [] RNA_5to3 = null;
	byte [] RNA_3to5 = null;
	int start;
	int stop;
	int length;
	
	public final static byte BASE_A = 0;
	public final static byte BASE_U = 1;
	public final static byte BASE_G = 2;
	public final static byte BASE_C = 3;
	
	//splice related
	private boolean [] forwardsStartLocations = null;
	private boolean [] forwardsBranchLocations = null;
	private boolean [] forwardsStopLocations = null;
	
	private boolean [] reverseStartLocations = null;
	private boolean [] reverseBranchLocations = null;
	private boolean [] reverseStopLocations = null;
	
	/**
	 * Assumes start < stop
	 * Assumes stop is an inclusive index
	 * @param dna
	 * @param start
	 * @param stop
	 */
	public RNA_Sequence (SequenceNucleotide sequenceFile, NucleotideSequence dna, int start, int stop) {
		this.sequenceFile = sequenceFile;
		DNA = dna;
		this.start = start;
		this.stop = stop;
		this.length = stop - start + 1;
		RNA_5to3 = new byte[length];
		RNA_3to5 = new byte[length];
		
		String DNAsequence = dna.getSequence();
		//convert the DNA to our forwards RNA strand
		for (int i = start; i < stop; i++) {
			RNA_5to3[i - start] = DNAtoRNA(DNAsequence.charAt(i));
		}
	
		//convert the RNA to the compliment
		for (int i = 0; i < length; i++) {
			//note that I'm putting it in reverse order
			RNA_3to5[i] = getRNACompliment(RNA_5to3[length - i - 1]);
		}
		
		//find our splices
		forwardsStartLocations = new boolean[length];
		for (int i = 0; i < length; i++) {forwardsStartLocations[i] = false;}
		forwardsBranchLocations = new boolean[length];
		for (int i = 0; i < length; i++) {forwardsBranchLocations[i] = false;}
		forwardsStopLocations = new boolean[length];
		for (int i = 0; i < length; i++) {forwardsStopLocations[i] = false;}
		
		reverseStartLocations = new boolean[length];
		for (int i = 0; i < length; i++) {reverseStartLocations[i] = false;}
		reverseBranchLocations = new boolean[length];
		for (int i = 0; i < length; i++) {reverseBranchLocations[i] = false;}
		reverseStopLocations = new boolean[length];
		for (int i = 0; i < length; i++) {reverseStopLocations[i] = false;}
		
		findSpliceLocations();
	}

	public SequenceNucleotide getSequenceFile() {
		return sequenceFile;
	}

	public static final double beliveableStartProbability = 0.01;
	public static final double beliveableBranchProbability = 0.00;
	public static final double beliveableStopProbability = 0.08;
	public static final double pyrimidineRichThresholdLevel = 0.725;
	public static final int pyrRichLength = 15;
	
	/**
	 *  
	 * @param dna assumes DNA is upper case
	 * @return
	 */
	public static byte DNAtoRNA(char dna) {
		if (dna == 'A') return BASE_A;
		if (dna == 'T') return BASE_U;
		if (dna == 'G') return BASE_G;
		return BASE_C;
	}
	
	/**
	 * A method for printing purposes
	 * @param rna
	 * @return
	 */
	public static char getRNAChar(byte rna) {
		if (rna == BASE_A) return 'A';
		if (rna == BASE_U) return 'U';
		if (rna == BASE_G) return 'G';
		return 'C';
	}
	
	public static byte getRNACompliment(byte rna) {
		if (rna == BASE_A) return BASE_U;
		if (rna == BASE_U) return BASE_A;
		if (rna == BASE_G) return BASE_C;
		return BASE_G;
	}
	
	public static boolean isPyrimidine(byte rna) {
		if (rna == BASE_C) return true;
		if (rna == BASE_U) return true;
		return false;
	}
	
	public void checkCD44Sites() {
//		int [] startZeroIndicies = {
//				501,
//				37871,
//				41538,
//				48031,
//				51196,
//				59377,
//				62326,
//				62918,
//				65771,
//				67374,
//				69337,
//				71185,
//				72580,
//				76045,
//				80518,
//				82863
//			};
//		int difference = -1;
//		int foundTotal = 0;
//		for (int i = 0; i < startZeroIndicies.length; i++) {
//			int foundSpot = startZeroIndicies[i] + difference;
//			if (forwardsStartLocations[foundSpot]) {
//				foundTotal++;
//			}
//			for (int j = foundSpot; j < foundSpot + 20; j++) {
//				System.out.print(getRNAChar(RNA_5to3[j]));
//			}
//			U.p();
//		}
//		U.p((double) foundTotal / startZeroIndicies.length);
		
		int [] stopZeroIndicies = {
				37705,
				41404,
				47962,
				50965,
				59251,
				62212,
				62801,
				65642,
				67242,
				69235,
				71095,
				72376,
				75982,
				80446,
				82784,
				90259
		};
		int difference = -2;
		int foundTotal = 0;
		for (int i = 0; i < stopZeroIndicies.length; i++) {
			int foundSpot = stopZeroIndicies[i] + difference;
			if (forwardsStopLocations[foundSpot]) {
				foundTotal++;
			}
			for (int j = foundSpot - 20; j <= foundSpot; j++) {
				System.out.print(getRNAChar(RNA_5to3[j]));
			}
			U.p();
		}
		U.p((double) foundTotal / stopZeroIndicies.length);

	}
	
	//http://www.ncbi.nlm.nih.gov/bookshelf/br.fcgi?book=mcb&part=A2868&rendertype=figure&id=A2883
	private void findSpliceLocations() {
		findSpliceSites(RNA_5to3, forwardsStartLocations, forwardsBranchLocations, forwardsStopLocations);
		findSpliceSites(RNA_3to5, reverseStartLocations, reverseBranchLocations, reverseStopLocations);
	}

	private void findSpliceSites(byte [] code, boolean [] startSites, boolean [] branchSites, boolean [] stopSites ) {
		int begin, end;
		double probability;
		
		//find start points
		begin = 3;
		end = length - 5;
		for (int i = begin; i < end; i++) {
			if (code[i] == BASE_G && code[i + 1] == BASE_U) {
				probability = 1.0;
				if (!(code[i -3] == BASE_A || code[i -3] == BASE_C )) probability *= 0.30;
				if (probability < beliveableStartProbability) continue;
				if (!(code[i -2] == BASE_A)) probability *= 0.40;
				if (probability < beliveableStartProbability) continue;
				if (!(code[i -1] == BASE_G)) probability *= 0.20;
				if (probability < beliveableStartProbability) continue;
				
				if (!(code[i +2] == BASE_A || code[i +2] == BASE_G )) probability *= 0.05;
				if (probability < beliveableStartProbability) continue;
				if (!(code[i +3] == BASE_A)) probability *= 0.30;
				if (probability < beliveableStartProbability) continue;
				if (!(code[i +4] == BASE_G)) probability *= 0.20;
				if (probability < beliveableStartProbability) continue;
				if (!(code[i +5] == BASE_U)) probability *= 0.55;
				if (probability < beliveableStartProbability) continue;
				
				startSites[i] = true;
			}
		}
		
		//find branch points
		begin = 3;
		end = length - 1;
		for (int i = begin; i < end; i++) {
			if (code[i] == BASE_A) {
				probability = 1.0;
				if (!(code[i -3] == BASE_C)) probability *= 0.20;
				if (probability < beliveableBranchProbability) continue;
				if (!(code[i -2] == BASE_U)) probability *= 0.10;
				if (probability < beliveableBranchProbability) continue;
				if (!(code[i -1] == BASE_A || code[i -1] == BASE_G )) probability *= 0.20;
				if (probability < beliveableBranchProbability) continue;
				
				if (!(code[i +1] == BASE_C || code[i +1] == BASE_U )) probability *= 0.20;
				if (probability < beliveableBranchProbability) continue;
				
				branchSites[i] = true;
				
			}
		}
		
		
		//find end points
		begin = 3;
		end = length - 1;
		for (int i = begin; i < end; i++) {
			if (code[i] == BASE_G && code[i - 1] == BASE_A) {
				probability = 1.0;
				if (!(code[i -2] == BASE_C)) probability *= 0.20;
				if (probability < beliveableStopProbability) continue;
				
				if (!(code[i +1] == BASE_G)) probability *= 0.40;
				if (probability < beliveableStopProbability) continue;
				
				stopSites[i] = true;
				
			}
		}
		
		//find pyrimidine rich areas upstream from stops
		//eliminate stops without these areas
		int minimumPyrAmount = (int) (pyrRichLength * pyrimidineRichThresholdLevel);
		for (int i = 0; i < length; i++) {
			if (stopSites[i]) {
				//keep track of our pyrimidines within a given window
				int pyrimidineTotal = 0;
				
				//the end of our window
				end = i - 50;
				if (end < 0) end = 0;
				
				//keep track of how far in window we have moved
				int windowDistance = 0;
				
				//did we find a rich area?
				boolean pyrAreaFound = false;
				
				//j moves upstream the length of our window
				int j = i;
				
				//moving backwards from stop point
				for (; j >= end; j--) {
					if (isPyrimidine(code[j])) {
						pyrimidineTotal++;
						if (pyrimidineTotal >= minimumPyrAmount) {
							pyrAreaFound = true;
//							U.p();
//							for (int k = j; k < j + pyrRichLength; k++) {
//								System.out.print(getRNAChar(code[k]));
//							}
							break;
						}
					}
					if (windowDistance >= pyrRichLength) {
						//subtract pyr if our pyrRichLength window has passed it
						if (isPyrimidine(code[j + pyrRichLength - 1])) {
							pyrimidineTotal--;
						}
					}
					
					windowDistance++;
				}
				
				//continue looking for branch points because we need one upstream from
				//the pyrimidine rich area
				boolean keepStop = false;
				if (pyrAreaFound) {
					for (; j >= end; j--) {
						if (branchSites[j]) {
							keepStop = true;
							break;
						}
					}
				}
				
				if (!keepStop) {
					stopSites[i] = false;
				}
			}
		}

		
		//eliminate lonely branch points
		for (int i = 0; i < length; i++) {
			if (branchSites[i]) {
				if (i > length - 51) {
					branchSites[i] = false;
					break;
				}
				boolean keepBranch = false;
				for (int j = i + 20; j <= i + 50; j++) {
					if (stopSites[j]) {
						keepBranch = true;
						break;
					}
				}
				if (!keepBranch) branchSites[i] = false;
			}
		}
		
		
		
	}

	public boolean[] getForwardsStartLocations() {
		return forwardsStartLocations;
	}

	public boolean[] getForwardsStopLocations() {
		return forwardsStopLocations;
	}

	public boolean[] getReverseStartLocations() {
		return reverseStartLocations;
	}

	public boolean[] getReverseStopLocations() {
		return reverseStopLocations;
	}

	public byte[] getRNA_3to5() {
		return RNA_3to5;
	}

	public byte[] getRNA_5to3() {
		return RNA_5to3;
	}

	public int getStart() {
		return start;
	}

	public int getStop() {
		return stop;
	}

	/**
	 * A wee function to see what kind of results we are getting.
	 */
	public void printStats() {
		U.p("sequence length: " + length);
		
		int startTotal = 0;
		for (int i = 0; i < length; i++) {
			if (forwardsStartLocations[i]) startTotal++;
		}
		U.p("Forwards start points: " + startTotal);
		
		int stopTotal = 0;
		for (int i = 0; i < length; i++) {
			if (forwardsStopLocations[i]) stopTotal++;
		}
		U.p("Forwards stop points: " + stopTotal);
		
		int prevLocation = 0;
		int count = 0;
		final int distance = 60;
		for (int i = 0; i < length; i++) {
			if (forwardsStopLocations[i]) {
				if (i - prevLocation < distance) {
					count++;
				}
				prevLocation = i;
			}
		}
		U.p("number of places where stops are less than " + distance + " apart: " + count);
		
		//printing a sample intron
		U.p("Example intron:");
		int intronStart = 0;
		int intronStop = 0;
		for (int i = 0; i < length; i++) {
			if (forwardsStartLocations[i]) {
				intronStart = i;
				break;
			}
		}
		for (int i = intronStart; i < length; i++) {
			if (forwardsStopLocations[i]) {
				intronStop = i;
				break;
			}
		}
		for (int i = intronStart; i <= intronStop; i++) {
			System.out.print(getRNAChar(RNA_5to3[i]));
		}
		U.p();
	}

}
