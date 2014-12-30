package Experimental;

import java.util.Enumeration;
import java.util.Hashtable;

import Navigator.Match;
import Navigator.MatchesToProtein;

/**
 * Many proteins contain similar sequences -- especially if they
 * derive from the same gene.  This will show how much overlap is between the two.
 * 
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class ProteinOverlap implements Comparable<ProteinOverlap> {
	
	MatchesToProtein proteinA;
	MatchesToProtein proteinB;
	
	int onlyACount = 0;
	int overlapCount = 0;
	int onlyBCount = 0;
	
	public ProteinOverlap(MatchesToProtein proteinA, MatchesToProtein proteinB) {
		this.proteinA = proteinA;
		this.proteinB = proteinB;
		
		Hashtable<String, Match> peptidesA = proteinA.getBestPeptides();
		Hashtable<String, Match> peptidesB = proteinB.getBestPeptides();
		
		Enumeration<String> e = peptidesA.keys();
		while (e.hasMoreElements()) {
			String peptide = e.nextElement();
			if (peptidesB.get(peptide) == null) {
				onlyACount++;
			} else {
				overlapCount++;
			}
		}
		
		onlyBCount = peptidesB.size() - overlapCount;
	}

	public int compareTo(ProteinOverlap other) {
		if (getOverlapCount() > other.getOverlapCount()) return -1;
		if (getOverlapCount() < other.getOverlapCount()) return  1;
		return 0;
	}

	public MatchesToProtein getProteinA() {
		return proteinA;
	}

	public MatchesToProtein getProteinB() {
		return proteinB;
	}

	public int getOnlyACount() {
		return onlyACount;
	}

	public int getOverlapCount() {
		return overlapCount;
	}

	public int getOnlyBCount() {
		return onlyBCount;
	}

}
