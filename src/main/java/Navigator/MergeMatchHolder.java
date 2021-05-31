package Navigator;

import java.util.ArrayList;
import java.util.Hashtable;

/**
 * Allows us to hold the best matches for all samples
 * <p>
 * Copyright 2013, Brian Risk
 *
 * @author Brian Risk
 */
public class MergeMatchHolder implements Comparable<MergeMatchHolder> {

    private String key = null;
    private double topScore = 0;
    private boolean isModified = true;
    private String keyName;
    private double scoreTotal = 0;

    Hashtable<String, MatchRow> matches = new Hashtable<String, MatchRow>();

    public MergeMatchHolder(String keyName) {
        this.keyName = keyName;
    }

    public void put(String key, MatchRow match) {
        matches.put(key, match);
        if (match.getScore() > topScore) topScore = match.getScore();
        scoreTotal += match.getScore();
        if (key == null) key = match.getString(keyName);
        if (!match.getBoolean("isModified")) isModified = false;
    }

    public MatchRow get(String key) {
        return matches.get(key);
    }

    public MatchRow get() {
        if (matches.size() == 0) return null;
        return matches.elements().nextElement();
    }

    public ArrayList<MatchRow> getMatches() {
        return new ArrayList<MatchRow>(matches.values());
    }

    public String getPeptideSequence() {
        return key;
    }


    public double getTopScore() {
        return topScore;
    }

    public double getScoreTotal() {
        return scoreTotal;
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
//		if (scoreTotal > other.getScoreTotal()) return -1;
//		if (scoreTotal < other.getScoreTotal()) return 1;
        return 0;
    }

}
