package Navigator;

import java.util.ArrayList;

/**
 * Copyright 2012, Brian Risk
 * Released under the Netscape Public License
 * 
 * @author Brian Risk
 *
 */
public class Region {
	
	String sequenceName;
	int start;
	int stop;
	double score;
	
	ArrayList<Match> matches = new ArrayList<Match>();
	
	
	public Region(String sequenceName, int start, int stop, double score) {
		super();
		this.sequenceName = sequenceName;
		this.start = start;
		this.stop = stop;
		this.score = score;
	}
	
	public void addMatch(Match match) {
		matches.add(match);
	}


	public String getSequenceName() {
		return sequenceName;
	}


	public void setSequenceName(String sequenceName) {
		this.sequenceName = sequenceName;
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



}
