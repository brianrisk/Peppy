package Navigator;

import java.util.ArrayList;

import Peppy.ModificationEntry;

/**
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class Modification implements Comparable<Modification> {
	
	double mass;
	Sample sample;
	
	double score = 0;
	
	ArrayList<Match> matches = new ArrayList<Match>();
	
	ArrayList<ModificationEntry> entries = new ArrayList<ModificationEntry>();
	
	public Modification(double mass, Sample sample) {
		this.mass = mass;
		this.sample = sample;
	}
	
	public Sample getSample() {
		return sample;
	}

	public void addMatch(Match match) {
		matches.add(match);
	}
	
	public void addToScore(double amount) {
		score += amount;
	}
	
	public  ArrayList<Match> getMatches() {
		return matches;
	}
	
	public void addEntry(ModificationEntry entry) {
		entries.add(entry);
	}
	
	public int getMatchesSize() {
		return matches.size();
	}
	
	public double getMass() {
		return mass;
	}

	
	public double getScore() {
		return score;
	}

	public int compareTo(Modification other) {
//		if (getMatchCount() > other.getMatchCount()) return -1;
//		if (getMatchCount() < other.getMatchCount()) return  1;
		if (getMass() > other.getMass()) return 1;
		if (getMass() < other.getMass()) return -1;
		return 0;
	}
	
	public ArrayList<ModificationEntry> getEntries() {
		return entries;
	}

}
