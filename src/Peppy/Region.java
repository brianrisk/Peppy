package Peppy;

import java.util.ArrayList;

/**
 * A way of grouping matches by regions and scoring those regions
 * 
 * initially this is just for DNA, but it could be for proteins
 * @author Brian Risk
 *
 */
public class Region implements Comparable<Region> {
	
	private ArrayList<Match> matches = new ArrayList<Match>();
	private Sequence sequence;
	private int startLocation;
	private int stopLocation;
	private int maxLength;
	private double score = 0;
	private int coverage = 0;
	private double pValue;
	private boolean flag = false;
	private boolean isForward;
	
	public Region(int startLocation, int maxLength, Sequence sequence, boolean isForward) {
		this.startLocation = startLocation;
		this.maxLength = maxLength;
		stopLocation = startLocation + maxLength;
		this.sequence = sequence;
		this.isForward = isForward;
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
//			if (match.rank == 1  && match.rankCount == 1) score += 1;
			if (match.rank == 1) score += (1.0 / match.rankCount);
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

	public int getMaxLength() {
		return maxLength;
	}

	public Sequence getSequence() {
		return sequence;
	}

	public double getPValue() {
		return pValue;
	}

	public double getScore() {
		return score;
	}
	
	public double calculateScore() {
//		score = 0;
//		for (Match match: matches) {
//			if (match.rank == 1) score += 1.0 / match.rankCount;
//		}
//		score = Math.log(score);
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

	public void setPValue(double pValue) {
		this.pValue = pValue;
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

	public boolean isForward() {
		return isForward;
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
		return matches.size();
	}
	
	public int getTrueStart() {
		int out = Integer.MAX_VALUE;
		for (Match match: matches) {
			if (match.getPeptide().getStartIndex() < out) out = match.getPeptide().getStartIndex();
		}
		return out;
	}
	
	public int getTrueStop() {
		int out = 0;
		for (Match match: matches) {
			if (match.getPeptide().getStopIndex() > out) out = match.getPeptide().getStartIndex();
		}
		return out;
	}
	
	/**
	 * returns the number of PSMs where rank = 1 and rankCount = 1
	 * @return
	 */
	public int getUniqueCount() {
		int out = 0;
		for (Match match: matches) {
			if (match.rank == 1 && match.rankCount == 1) out++;
		}
		return out;
	}
}
