package Navigator;

import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;


/**
 * Contains a hashtable where the elements are
 * ArrayList of matches.
 * 
 * Copyright 2012, Brian Risk
 * Released under the Netscape Public License
 * 
 * @author Brian Risk
 *
 */
public class MatchTable { 
	
	/* where we store all the matches */
	private Hashtable<String, ArrayList<Match>> matches = new Hashtable<String, ArrayList<Match>>();
	
	/*
	 * This is a flag that sets if we are to clear out the existing ArrayList of matches if
	 * we are adding a match of greater score.
	 */
	private boolean keepOnlyTheBest;
	
	/* in some cases we will want to keep only one match per key */
	
	
	public MatchTable(boolean keepOnlyTheBest) {
		this.keepOnlyTheBest = keepOnlyTheBest;
	}
	
	
	/**
	 * Handles all of the headache of adding to a hashtable of array lists.
	 * @param key
	 * @param match
	 */
	public void put(String key, Match match) {
		ArrayList<Match> existingMatches = matches.get(key);
		if (existingMatches == null) {
			existingMatches = new ArrayList<Match>();
			existingMatches.add(match);
			matches.put(key, existingMatches);
		} else {

			
			/*
			 * if we are keeping only the best then we must see what the score
			 * of the existing matches are (and they should all be the same if
			 * there are more than one).  
			 */
			if (keepOnlyTheBest) {
				
				/* If the new match has a higher score then we create a new ArrayList and usurp the existing with it. */
				if (match.getScore() > existingMatches.get(0).getScore()) {
					ArrayList<Match> matchArray = new ArrayList<Match>();
					matchArray.add(match);
					matches.put(key, matchArray);
				}
				
				/* If the new score equals the existing score, we add it to the list */
				if (match.getScore() == existingMatches.get(0).getScore()) {
					existingMatches.add(match);
				}
				
			} else {
				existingMatches.add(match);
			}
		}

	}
	
	
	/**
	 * 
	 * @return
	 */
	public ArrayList<Match> getArrayList() {
		/* determining the size of the output */
		int size = 0;
		Enumeration<ArrayList<Match>> e = matches.elements();
		while (e.hasMoreElements()) {
			size += e.nextElement().size();
		}
		
		/* create the output */
		ArrayList<Match> out = new ArrayList<Match>(size);
		e = matches.elements();
		while (e.hasMoreElements()) {
			out.addAll(e.nextElement());
		}
		
		return out;
	}
	
	public Hashtable<String, ArrayList<Match>> getHashtable() {
		return matches;
	}
	
	public ArrayList<Match> get(String key) {
		return matches.get(key);
	}

	
}
