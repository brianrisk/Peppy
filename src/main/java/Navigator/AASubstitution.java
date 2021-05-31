package Navigator;

import Peppy.AminoAcids;

import java.util.ArrayList;
import java.util.Collections;

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
