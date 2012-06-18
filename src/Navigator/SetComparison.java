package Navigator;

import java.util.Enumeration;
import java.util.Hashtable;

/**
 * Copyright 2012, Brian Risk
 * Released under the Netscape Public License
 * 
 * @author Brian Risk
 *
 */
public class SetComparison implements Comparable<SetComparison>{
	
	private String name;
	
	private Hashtable<String, MatchesTo> matchesTos = new Hashtable<String, MatchesTo>();
	
	/* the comparison score is the largest individual / average ratio */
	private double score = 0;
	
	private double total = 0;
	
	private double maxIndividualScore = 0;
	
	private String dominantSetName = "";
	
	private MatchesTo dominantSet;
	
	private int maxMatchCount = 0;
	
	private int minMatchCount = Integer.MAX_VALUE;
	
	
	
	public SetComparison(String name) {
		this.name = name;
	}
	
	public void addSet(String sampleName, MatchesTo matchesTo) {
		matchesTos.put(sampleName, matchesTo);
		total += matchesTo.getScore();
		
		if (matchesTo.getMatchesSize() > maxMatchCount) maxMatchCount = matchesTo.getMatchesSize();
		
		if (matchesTo.getMatchesSize() < minMatchCount) minMatchCount = matchesTo.getMatchesSize();
		
		
		/* maybe this new score is better */
		if (maxIndividualScore < matchesTo.getScore()) {
			maxIndividualScore = matchesTo.getScore();
			dominantSetName = sampleName;
			dominantSet = matchesTo;
		}	
		
		/* recalculate our score */
		score = getRelativeScore(maxIndividualScore);		
		
	}


	public double getRelativeScore(String sampleName) {
		double sampleScore = matchesTos.get(sampleName).getScore();
		return getRelativeScore(sampleScore);
	}
	
	public double getRelativeScore(double score) {
		if (score == total) return score;
		double averageMinusTop =  (total - score) /(matchesTos.size() - 1);
		if (averageMinusTop == 0) {
			return Double.MAX_VALUE;
		} else {
			return score / averageMinusTop;
		}
	}


	public int compareTo(SetComparison other) {
		if (getScore() > other.getScore()) return -1;
		if (getScore() < other.getScore()) return  1;
		return 0;
	}


	public int getMaxMatchCount() {
		return maxMatchCount;
	}

	public int getMinMatchCount() {
		return minMatchCount;
	}

	public Hashtable<String, MatchesTo> getSets() {
		return matchesTos;
	}


	public double getScore() {
		return score;
	}

	public String getName() {
		return name;
	}

	public String getDominantSampleName() {
		return dominantSetName;
	}

	public MatchesTo getDominantSet() {
		return dominantSet;
	}
	
	public int getMatchTotal() {
		int out = 0;
		Enumeration<MatchesTo> e = matchesTos.elements();
		while (e.hasMoreElements()) {
			MatchesTo matches = e.nextElement();
			out += matches.getMatchesSize();
		}
		return out;
	}
	

}
