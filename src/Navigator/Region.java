package Navigator;

import java.util.ArrayList;

/**
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class Region {
	
	String name = "";
	String sequence = "";
	String description = "";
	int start = -1;
	int stop = -1;
	boolean isForwards = true;
	double score = 0;
	
	ArrayList<Match> matches = new ArrayList<Match>();
	
	
	public Region(String name, int start, int stop, boolean isForwards) {
		super();
		this.name = name;
		this.start = start;
		this.stop = stop;
	}
	
	public Region(){}
	
	public void addMatch(Match match) {
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



}
