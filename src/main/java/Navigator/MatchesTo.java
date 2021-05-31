package Navigator;

import java.util.ArrayList;

/**
 * The MatchesTo superset is a collection of matches.
 * These could be all matches to a particular peptide or
 * all for a given region or protein or so on and so on.
 * <p>
 * The convenience of this superclass is that it treats the
 * top score for all children matches as the score for the set.
 * <p>
 * This allows sorting of match sets.
 * <p>
 * Copyright 2013, Brian Risk
 *
 * @author Brian Risk
 */
public abstract class MatchesTo implements Comparable<MatchesTo> {


    private double score = 0;

    private String name;

    private ArrayList<MatchRow> matches = new ArrayList<MatchRow>();


    public MatchesTo(String name) {
        this.name = name;
    }

    public void addMatch(MatchRow match) {
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

    public ArrayList<MatchRow> getMatches() {
        return matches;
    }

    public int getMatchesSize() {
        return matches.size();
    }

    public int compareTo(MatchesTo other) {
        if (score > other.getScore()) return -1;
        if (score < other.getScore()) return 1;
        return 0;
    }


}
