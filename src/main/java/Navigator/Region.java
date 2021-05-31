package Navigator;

import java.util.ArrayList;

/**
 * Copyright 2013, Brian Risk
 *
 * @author Brian Risk
 */
public class Region {

    String name = "";
    String sequence = "";
    String description = "";
    int start = -1;
    int stop = -1;
    boolean isForwards = true;
    double score = 0;

    ArrayList<MatchRow> matches = new ArrayList<MatchRow>();

    public Region() {
    }

    public void addMatch(MatchRow match) {
        matches.add(match);
    }


    public String getName() {
        return name;
    }


    public void setName(String sequenceName) {
        this.name = sequenceName;
    }


    public int getStart() {
        return start;
    }


    public void setStart(int start) {
        this.start = start;
    }


    public int getStop() {
        return stop;
    }


    public void setStop(int stop) {
        this.stop = stop;
    }

    /**
     * gets the width that this region covers
     *
     * @return
     */
    public int getCoverage() {
        return stop - start;
    }


    public double getScore() {
        return score;
    }


    public void setScore(double score) {
        this.score = score;
    }

    public boolean isForwards() {
        return isForwards;
    }

    public void setForwards(boolean isForwards) {
        this.isForwards = isForwards;
    }

    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    /**
     * This method build for the code that will turn a GTF sequence
     * into a six-frame translation with some padding around the transcript
     * boundaries.
     *
     * @param nucleotides
     */
    public void addPadding(int nucleotides) {
        start -= nucleotides;
        if (start < 0) start = 0;
        stop += nucleotides;
    }

    public boolean intersects(Region region) {
        if (region.getSequence().equals(sequence)) {
            if (region.getStart() >= start && region.getStart() < stop) return true;
            if (region.getStop() >= start && region.getStop() < stop) return true;
        }
        return false;
    }

    public boolean addRegion(Region region) {
        if (intersects(region)) {
            if (region.start < start) start = region.start;
            if (region.stop > stop) stop = region.stop;
            return true;
        } else {
            return false;
        }
    }


}
