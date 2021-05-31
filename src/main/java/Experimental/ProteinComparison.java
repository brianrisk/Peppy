package Experimental;

import Navigator.MatchesTo;
import Navigator.MatchesToProtein;

import java.util.Enumeration;

/**
 * Copyright 2013, Brian Risk
 *
 * @author Brian Risk
 */
public class ProteinComparison extends SetComparison {


    private int maxUniquePeptideCount = 0;

    private double maxUniquePeptideScore = 0;

    private boolean hasUnique = false;


    public ProteinComparison(String proteinName) {
        super(proteinName);
    }

    public void addSet(String sampleName, MatchesToProtein matchesToProtein) {
        super.addSet(sampleName, matchesToProtein);

        if (matchesToProtein.getUniquePeptideCount() > maxUniquePeptideCount)
            maxUniquePeptideCount = matchesToProtein.getUniquePeptideCount();

        if (matchesToProtein.getMaxUniquePeptideScore() > maxUniquePeptideScore)
            maxUniquePeptideScore = matchesToProtein.getMaxUniquePeptideScore();

        if (matchesToProtein.hasUniqueMatch()) hasUnique = true;

    }

    public boolean hasUnique() {
        return hasUnique;
    }


    public int getMaxUniquePeptideCount() {
        return maxUniquePeptideCount;
    }

    public double getMaxUniquePeptideScore() {
        return maxUniquePeptideScore;
    }

    public boolean isFullyOverlapped() {
        Enumeration<MatchesTo> e = getSets().elements();
        while (e.hasMoreElements()) {
            if (!((MatchesToProtein) e.nextElement()).isEntirelyFoundInOtherProtein()) return false;
        }
        return true;
    }

}
