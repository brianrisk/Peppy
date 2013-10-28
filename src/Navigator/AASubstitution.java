package Navigator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

import Peppy.AminoAcids;
import Peppy.DigestionThread_DNA;

/**
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class AASubstitution implements Comparable<AASubstitution>{
	
	/* the former and the present amino acids */
	private char previous;
	private char present;
	private double difference;
	
	/**
	 * this constructor is used to make a blank object that can be used when
	 * searching through a sorted list.
	 * @param difference
	 */
	public AASubstitution(double difference) {
		this.difference = difference;
	}
	
	public AASubstitution(char previous, char present, double difference) {
		this.previous = previous;
		this.present = present;
		this.difference = difference;
	}
	

	
	public static ArrayList<AASubstitution> generateListOfAASubstitutions() {
		ArrayList<AASubstitution> out = new ArrayList<AASubstitution>(); 
		char [] acids = AminoAcids.getAcids();

		/* starting with 1 to skip the stop */
		double difference;
		for (int i = 1; i < acids.length; i++) {
			for (int j = 1; j < acids.length; j++) {
				difference = AminoAcids.getWeightMono(acids[i]) - AminoAcids.getWeightMono(acids[j]);
				if (Math.abs(difference) > 0.01) {
					out.add(new AASubstitution(acids[i], acids[j], -difference));
				}
			}
		}
		
		Collections.sort(out);
		
		return out;

	}
	
	
	public static ArrayList<AASubstitution> generateListOfSingleNucleotideAASubstitutions() {
		Hashtable<String, AASubstitution> substitutionHash = new Hashtable<String, AASubstitution>();
		char [] nucleotides = {'A', 'T', 'G', 'C'};
		char [] codon = new char[3];
		char [] codon2 = new char[3];
		char acid, acid2;
		double difference;
		for (int indexA = 0; indexA < 4; indexA++) {
			for (int indexB = 0; indexB < 4; indexB++) {
				for (int indexC = 0; indexC < 4; indexC++) {
					codon[0] = nucleotides[indexA];
					codon[1] = nucleotides[indexB];
					codon[2] = nucleotides[indexC];
					acid = AminoAcids.aminoAcidList[DigestionThread_DNA.indexForCodonArray(codon)];
					
					for (int codonIndex = 0; codonIndex < 3; codonIndex++) {
						for (int snpIndex = 0; snpIndex < 4; snpIndex++) {
							codon2[0] = nucleotides[indexA];
							codon2[1] = nucleotides[indexB];
							codon2[2] = nucleotides[indexC];
							/* the SNP */
							codon2[codonIndex] = nucleotides[snpIndex];
							acid2 = AminoAcids.aminoAcidList[DigestionThread_DNA.indexForCodonArray(codon2)];
							
							if (acid != acid2) {
								difference = AminoAcids.getWeightMono(acid) - AminoAcids.getWeightMono(acid2);
								if (Math.abs(difference) > 0.01) {
									AASubstitution aasub = new AASubstitution(acid, acid2, -difference);
									String hash = acid + "-" + acid2;
									substitutionHash.put(hash, aasub);
								}
								
							}
						}
					}
				}
			}
		}
		
		return new ArrayList<AASubstitution>(substitutionHash.values());
	}



	public char getPrevious() {
		return previous;
	}



	public char getPresent() {
		return present;
	}



	public double getDifference() {
		return difference;
	}



	public int compareTo(AASubstitution o) {
		if (getDifference() > o.getDifference()) return 1;
		if (getDifference() < o.getDifference()) return -1;
		return 0;
	}

}
