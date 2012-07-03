package Peppy;

import java.util.ArrayList;
import java.util.Collections;

/**
 * We want to hold collections of matches for various things.  All matches for a given peptide,
 * or spectrum or region or protein... This abstract class defines common functionality for all
 * of these.
 * 
 * Copyright 2012, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public abstract class Matches implements Comparable<Matches> {
	
	ArrayList<Match> matches = new ArrayList<Match>();
	
	/* the top score of all of our matches */
	double score = Properties.minimumScore;
		
	
	private static int labelTracker = 0;
	public static final int KEEP_ONLY_BEST_MATCHES = labelTracker++;
	public static final int KEEP_MATCHES_AT_MINIMUM_SCORE = labelTracker++;
	
	private int whatToKeep = KEEP_ONLY_BEST_MATCHES;
	
	/* only keep decoy hits if they trump our existing score 
	 * For FDR we only need to know when the decoy trumps the target hits */
	private boolean ignoreLesserDecoys = true;
	
	/**
	 * Add a match only if it's score is greater than or equal to the reigning score.
	 * If it is greater than, then the existing matches are cleared out.
	 * 
	 * NOTE: don't add if it is a duplicate match??
	 * 
	 * @param match
	 */
	public void addMatch(Match match) {
		if (ignoreLesserDecoys) {
			if (match.getPeptide().isDecoy()) {
				if (match.getScore() <= score) {
					return;
				} 
			}
		}
		
		if (whatToKeep == KEEP_ONLY_BEST_MATCHES) {
			if (match.getScore() > score) {
				matches.clear();
				matches.add(match);
				score = match.getScore();
			} else {
				if (match.getScore() == score) {
					matches.add(match);
				}
			}
			return;
		}
		

		if (whatToKeep == KEEP_MATCHES_AT_MINIMUM_SCORE) {
			if (match.getScore() >= Properties.minimumScore) {
				matches.add(match);
				if (match.getScore() > score) {
					score = match.getScore();
				}
			} 
		}
	}
	
	public void clearMatches() {
		matches.clear();
		score = Properties.minimumScore;
	}
	
	public ArrayList<Match> getMatches() {
		return matches;
	}

	public double getScore() {
		if (matches.size() > 0) {
			return score;
		} else {
			return 0;
		}
	}
	
	public static ArrayList<Match> getBestMatches(ArrayList<Match> matches) {
		ArrayList<Match> bestMatches = new ArrayList<Match>();
		Match.setSortParameter(Match.SORT_BY_PEPTIDE_THEN_SCORE);
		Collections.sort(matches);
		Peptide peptide = new Peptide("k");
		for (Match match: matches) {
			if (!peptide.equals(match.getPeptide())) {
				peptide = match.getPeptide();
				bestMatches.add(match);
			}
		}
		Match.setSortParameter(Match.SORT_BY_SCORE);
		Collections.sort(matches);
		Collections.sort(bestMatches);
		return bestMatches;
	}
	
	public static ArrayList<Match> getMatchesWithSpectrum(Spectrum spectrum, ArrayList<Match> theseMatches) {
		ArrayList<Match> out = new ArrayList<Match>();
		for (int i = 0; i < theseMatches.size(); i++) {
			Match match = theseMatches.get(i);
			if (match.getSpectrum().getId() == spectrum.getId()) {
				out.add(match);
			}
		}
		return out;
	}
	
	public static ArrayList<Match> getMatchesWithSequence(Sequence_DNA sequence_DNA, ArrayList<Match> theseMatches) {
		ArrayList<Match> out = new ArrayList<Match>();
		for (int i = 0; i < theseMatches.size(); i++) {
			Match match = theseMatches.get(i);
			if (match.getPeptide().getParentSequence().getId() == sequence_DNA.getId()) {
				out.add(match);
			}
		}
		return out;
	}

	public int getWhatToKeep() {
		return whatToKeep;
	}

	public void setWhatToKeep(int whatToKeep) {
		this.whatToKeep = whatToKeep;
	}
	
	
	public int compareTo(Matches o) {
		if (getScore() < o.getScore()) return 1;
		if (getScore() > o.getScore()) return -1;
		return 0;
	}

	public boolean ignoreLesserDecoys() {
		return ignoreLesserDecoys;
	}

	public void setIgnoreLesserDecoys(boolean ignoreLesserDecoys) {
		this.ignoreLesserDecoys = ignoreLesserDecoys;
	}

}
