package Peppy;

import java.util.ArrayList;
import java.util.Collections;

public class Matches {
	
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
		Match.setSortParameter(Match.SORT_BY_IMP_VALUE);
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

}
