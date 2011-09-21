package Peppy;

import java.util.ArrayList;

import Math.HasEValue;

/**
 * A way of grouping matches by regions and scoring those regions
 * 
 * initially this is just for DNA, but it could be for proteins
 * @author Brian Risk
 *
 */
public class Region implements HasEValue, Comparable<Region> {
	
	private ArrayList<Match> matches = new ArrayList<Match>();
	private Sequence sequence;
	private int startLocation;
	private int stopLocation;
	private int maxLength;
	private double score = 0;
	private int coverage = 0;
	private int numberOfMatches = 0;
	private double eValue;
	private boolean flag = false;
	
	public Region(int startLocation, int maxLength, Sequence sequence) {
		this.startLocation = startLocation;
		this.maxLength = maxLength;
		stopLocation = startLocation + maxLength;
		this.sequence = sequence;
	}
	
	/**
	 * Adds a match if it starts within the regions boundaries
	 * @param match
	 * @return returns true if match was successfully added
	 */
	public boolean addMatch(Match match) {
		if (match.getPeptide().getStartIndex() >= startLocation && match.getPeptide().getStartIndex() < stopLocation) {
			matches.add(match);
//			score += match.getPeptide().getLength() * match.getPeptide().getLength();
//			score -= Math.log10(match.getEValue());
			score += 1;
			numberOfMatches++;
			return true;
		} else {
			return false;
		}
	}

	public int getStopLocation() {
		return stopLocation;
	}

	public int getStartLocation() {
		return startLocation;
	}

	public Sequence getSequence() {
		return sequence;
	}

	public double getEValue() {
		return eValue;
	}

	public double getScore() {
		return score;
	}
	
	public double calculateScore() {
		score = Math.log(score);
		return score;
	}
	
	public int calculateCoverage() {
		boolean[] coverageArray = new boolean[maxLength];
		int i, start, stop;
		for (Match match: matches) {
			start = match.getPeptide().getStartIndex() - startLocation;
			stop = match.getPeptide().getStopIndex() - startLocation;
			if (stop > maxLength) stop = maxLength;
			for (i = start; i < stop; i++) {
				coverageArray[i] = true;
			}
		}
		coverage = 0;
		for (i = 0; i < maxLength; i++) {
			if (coverageArray[i]) coverage++;
		}
		return coverage;
	}

	public void setEValue(double eValue) {
		this.eValue = eValue;
	}


	public int compareTo(Region o) {
		if (score < o.getScore()) return 1;
		if (score > o.getScore()) return -1;
		return 0;
	}

	public ArrayList<Match> getMatches() {
		return matches;
	}
	
	public boolean isOverlapping(Region other) {
		return (isWithin(other.getStartLocation()) || isWithin(other.getStopLocation()));
	}
	
	public boolean isWithin(int index) {
		return (index >= startLocation && index < stopLocation);
	}
	
	/* flagging */
	public boolean isFlagged() {
		return flag;
	}
	
	public boolean isUnFlagged() {
		return !flag;
	}

	public void flag() {
		flag = true;
	}
	
	public void unFlag() {
		flag = false;
	}
	
	public void clearRegion() {
		matches.clear();
	}
	
	public int getNumberOfMatches() {
		return numberOfMatches;
	}
}
