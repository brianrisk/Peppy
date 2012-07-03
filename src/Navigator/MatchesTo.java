package Navigator;

import java.util.ArrayList;

/**
 * A "set" is a collection of matches.  These could be all matches to a particular peptide or 
 * all for a given region or protein or so on and so on.
 * 
 * Copyright 2012, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public abstract class MatchesTo implements Comparable<MatchesTo>{

	
	private double score = 0;
	
	private String name;
	
	private ArrayList<Match> matches = new ArrayList<Match>();
	


	public MatchesTo(String name) {
		this.name = name;
	}
	
	public void addMatch(Match match) {
		matches.add(match);	
	}
	
	public void addToScore(double amount) {
		score += amount;
	}

	public double getScore() {
		return score;
	}

	public String getName() {
		return name;
	}

	public ArrayList<Match> getMatches() {
		return matches;
	}
	
	public int getMatchesSize() {
		return matches.size();
	}

	public int compareTo(MatchesTo other) {
		if (score > other.getScore()) return -1;
		if (score < other.getScore()) return  1;
		return 0;
	}
	




}
