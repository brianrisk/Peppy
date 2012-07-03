package Navigator;

import java.util.Hashtable;

/**
 * Allows us to hold the best matches for all samples
 * 
 * Copyright 2012, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class MergeMatchHolder implements Comparable<MergeMatchHolder> {
	
	private String key = null;
	private double topScore = 0;
	private boolean isModified = true;
	private String keyName;
	
	Hashtable<String, Match> matches = new Hashtable<String, Match>();
	
	public MergeMatchHolder(String keyName) {
		this.keyName = keyName;
	}
	
	public void put(String key, Match match) {
		matches.put(key, match);
		if (match.getScore() > topScore) topScore = match.getScore();
		if (key == null) key = match.getString(keyName);
		if (!match.getBoolean("isModified")) isModified = false;
	}
	
	public Match get(String key) {
		return matches.get(key);
	}
	
	public Match get() {
		if (matches.size() == 0) return null;
		return matches.elements().nextElement();
	}
	
	public String getPeptideSequence() {
		return key;
	}


	public double getTopScore() {
		return topScore;
	}
	
	public int size() {
		return matches.size();
	}
	
	public boolean isModified() {
		return isModified;
	}


	/* sorts first by size then score */
	public int compareTo(MergeMatchHolder other) {
		if (size() > other.size()) return -1;
		if (size() < other.size()) return 1;
		if (topScore > other.getTopScore()) return -1;
		if (topScore < other.getTopScore()) return 1;
		return 0;
	}

}
