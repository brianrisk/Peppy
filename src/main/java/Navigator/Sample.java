package Navigator;

import java.io.File;
import java.util.ArrayList;

/**
 * Copyright 2013, Brian Risk
 *
 * @author Brian Risk
 */
public class Sample {

    String name;

    private static int analysisDatabase = 0;
    public static final int CONTAMINANT_PROTEIN = analysisDatabase++;
    public static final int REFERENCE_PROTEIN = analysisDatabase++;
    public static final int SUBJECT_PROTEIN = analysisDatabase++;
    public static final int DISEASE_PROTEIN = analysisDatabase++;

    public static final int CONTAMINANT_GENOME = analysisDatabase++;
    public static final int REFERENCE_GENOME = analysisDatabase++;
    public static final int SUBJECT_GENOME = analysisDatabase++;
    public static final int DISEASE_GENOME = analysisDatabase++;

    private int databaseType;

    /* where we hold our best matches.
     * Sending it "true" meaning, keep only the best matches */
    MatchTable bestTable = new MatchTable(true);

    /* where we hold all of the matches we have loaded */
    ArrayList<MatchRow> matches = new ArrayList<MatchRow>();

    /* constructor */
    public Sample(String name, int databaseType) {
        this.name = name;
        this.databaseType = databaseType;
    }

    public void loadResults(File file) {
        ArrayList<MatchRow> newMatches = MatchRow.loadMatches(file);
        matches.addAll(newMatches);
        addMatches(newMatches);
    }

    public String getName() {
        return name;
    }

    public int getDatabaseType() {
        return databaseType;
    }

    public ArrayList<MatchRow> getMatches() {
        return matches;
    }

    /* add a set of matches to our MatchTable */
    public void addMatches(ArrayList<MatchRow> matches) {
        for (MatchRow match : matches) match.set("analysisDatabase", databaseType);

        /* is this from a protein or genome? */
        boolean isProtein = true;
        if (
                analysisDatabase == CONTAMINANT_GENOME ||
                        analysisDatabase == REFERENCE_GENOME ||
                        analysisDatabase == SUBJECT_GENOME ||
                        analysisDatabase == DISEASE_GENOME
                ) {
            isProtein = false;
        }

        /* label each of the matches what we found above */
        if (isProtein) {
            for (MatchRow match : matches) match.set("databaseType", "protein");
        } else {
            for (MatchRow match : matches) match.set("databaseType", "genome");
        }


        /* add all the matches */
        for (MatchRow match : matches) bestTable.put(match.getString("spectrumMD5"), match);
    }

    public void reduceMatches(double matchScoreCutoff) {
        int keepCount = 0;
        for (MatchRow match : matches) {
            if (match.getScore() >= matchScoreCutoff) keepCount++;
        }
        ArrayList<MatchRow> reducedMatches = new ArrayList<MatchRow>(keepCount);
        for (MatchRow match : matches) {
            if (match.getScore() >= matchScoreCutoff) reducedMatches.add(match);
        }
        matches = reducedMatches;
    }

}
