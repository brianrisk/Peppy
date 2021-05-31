

package Navigator;

import java.util.ArrayList;

/**
 * Copyright 2013, Brian Risk
 *
 * @author Brian Risk
 */
public class ModificationComparison implements Comparable<ModificationComparison> {

    String name;

    private ArrayList<Modification> modifications = new ArrayList<Modification>();

    /* the comparison score is the largest individual / average ratio */
    private double score = 0;

    private double totalScore = 0;

    private String dominantSampleName = "";

    private int maxMatchCount = 0;

    private int minMatchCount = Integer.MAX_VALUE;

    private double maxIndividualScore = 0;


    public ModificationComparison(String name) {
        this.name = name;
    }


    public void addModification(Modification modification) {
        modifications.add(modification);

        if (modification.getMatchesSize() > maxMatchCount) maxMatchCount = modification.getMatchesSize();

        if (modification.getMatchesSize() < minMatchCount) minMatchCount = modification.getMatchesSize();

        /* maybe this new score is better */
        if (maxIndividualScore < modification.getScore()) {
            maxIndividualScore = modification.getScore();
            dominantSampleName = modification.getSample().getName();
        }

        /* recalculate our score */
        totalScore += modification.getScore();
        score = getRelativeScore(maxIndividualScore);
    }


    public double getRelativeScore(double score) {
        if (score == totalScore) return score;
        double averageMinusTop = (totalScore - score) / (modifications.size() - 1);
        if (averageMinusTop == 0) {
            return Double.MAX_VALUE;
        } else {
            return score / averageMinusTop;
        }
    }


    public int compareTo(ModificationComparison other) {
        if (getScore() > other.getScore()) return -1;
        if (getScore() < other.getScore()) return 1;
        return 0;
    }


    public int getMaxMatchCount() {
        return maxMatchCount;
    }

    public int getMinMatchCount() {
        return minMatchCount;
    }


    public double getScore() {
        return score;
    }

    public String getName() {
        return name;
    }

    public String getDominantSampleName() {
        return dominantSampleName;
    }

    public double getMaxIndividualScore() {
        return maxIndividualScore;
    }


}
