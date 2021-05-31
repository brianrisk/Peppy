package Navigator;

import Math.MathFunctions;

/**
 * Used in "ModificationReport"
 * <p>
 * This object scores the probability that a modification
 * of a certain mass of a certain acid is random.  This
 * is computed by finding the frequency of the acids, the number
 * of times the modification took place and the number of times
 * it took place on this particular acid.
 *
 * @author Brian Risk
 */

public class ModificationType implements Comparable<ModificationType> {

    char acid;
    int mass;
    int tally = 0;
    double score;

    public ModificationType(char acid, int mass) {
        this.acid = acid;
        this.mass = mass;
    }


    public void calculateScore(int totalModOccurrences, double acidPercentage) {
        score = MathFunctions.approximateBinomialProbability(totalModOccurrences, tally, acidPercentage);
        score = -Math.log10(score);
        if (Double.isInfinite(score)) score = 500 + tally;
    }

    public void incrementTally() {
        tally++;
    }


    public int compareTo(ModificationType arg0) {
        if (score < arg0.getScore()) return 1;
        if (score > arg0.getScore()) return -1;
        return 0;
    }


    public char getAcid() {
        return acid;
    }


    public int getMass() {
        return mass;
    }


    public double getScore() {
        return score;
    }


    public int getTally() {
        return tally;
    }

}
