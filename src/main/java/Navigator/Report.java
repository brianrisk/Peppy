package Navigator;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

/**
 * An abstract class that objects that reports should inherit
 *
 * @author Brian Risk
 */
public abstract class Report {
    ArrayList<MatchRow> matches;


    /**
     * all reports are based on matches.  Thus, the constructor takes a set of matches
     *
     * @param matches
     */
    public Report(ArrayList<MatchRow> matches) {
        this.matches = matches;
    }


    /**
     * The file you return is the file or directory that this method creates
     *
     * @param reportDirectory
     * @return
     */
    public abstract File createReport(File reportDirectory) throws IOException;


    public void addMatches(ArrayList<MatchRow> matches) {
        this.matches.addAll(matches);
    }


    /**
     * A utility class to load in all matches from a report
     *
     * @param reportFile
     * @return
     */
    public static ArrayList<MatchRow> loadMatches(File reportFile) {
        return loadMatches(reportFile, -1);
    }


    /**
     * A utility class to load in matches and set a lowest score
     *
     * @param reportFile
     * @param scoreCutoff
     * @return
     */
    public static ArrayList<MatchRow> loadMatches(File reportFile, double scoreCutoff) {
        ArrayList<MatchRow> matches = MatchRow.loadMatches(reportFile);
        if (scoreCutoff > 0) {
            matches = filter(matches, scoreCutoff);
        }
        return matches;
    }


    private static ArrayList<MatchRow> filter(ArrayList<MatchRow> matches, double scoreCutoff) {
        ArrayList<MatchRow> out = new ArrayList<MatchRow>(matches.size());
        for (MatchRow match : matches) {
            if (match.getScore() < scoreCutoff) continue;
            out.add(match);
        }
        return out;
    }
}
