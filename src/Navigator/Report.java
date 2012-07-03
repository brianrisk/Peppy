package Navigator;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

/**
 * An abstract class that objects that create reports should inherit
 * 
 * @author Brian Risk
 *
 */
public abstract class Report {
	ArrayList<Match> matches;
	
	
	/**
	 * all reports are based on matches.  Thus, the constructor takes a set of matches
	 * @param matches
	 */
	public Report(ArrayList<Match> matches) {
		this.matches = matches;
	}
	
	
	
	/**
	 * The file you return is the file or directory that this method creates
	 * @param reportDirectory
	 * @return
	 */
	public abstract File createReport(File reportDirectory) throws IOException;
	
	
	
	
	public void addMatches(ArrayList<Match> matches) {
		this.matches.addAll(matches);
	}



	/**
	 * A utility class to load in all matches from a report
	 * @param reportFile
	 * @return
	 */
	public static ArrayList<Match> loadMatches(File reportFile) {
		return loadMatches(reportFile, -1);
	}

	
	/**
	 * A utility class to load in matches and set a lowest score
	 * @param reportFile
	 * @param scoreCutoff
	 * @return
	 */
	public static ArrayList<Match> loadMatches(File reportFile, double scoreCutoff) {
		ArrayList<Match> matches = Match.loadMatches(reportFile);
		if (scoreCutoff > 0) {
			matches = filter(matches, scoreCutoff);
		}
		return matches;
	}
	

	private static ArrayList<Match> filter(ArrayList<Match> matches, double scoreCutoff) {
		ArrayList<Match> out = new ArrayList<Match>(matches.size());
		for (Match match: matches) {
			if (match.getScore() < scoreCutoff) continue;
			out.add(match);
		}
		return out;
	}
}
