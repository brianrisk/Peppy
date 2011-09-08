package Reports;

import java.util.ArrayList;

import Peppy.Match;
import Peppy.Sequence_DNA;
import Peppy.Spectrum;

/**
 * These are common searches that one might perform on a body of matches
 * @author Brian Risk
 *
 */

public class CommonMatchSearches {
	
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
