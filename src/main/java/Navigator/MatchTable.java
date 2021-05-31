package Navigator;

import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;


/**
 * 
 * For any given peptide acid sequence, there can be a plethora of PSMs.
 * This data structure allows easy grouping by acid sequence.
 * 
 * Contains a hashtable where the elements are
 * ArrayList of matches.
 * 
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class MatchTable { 
	
	/* where we store all the matches */
	private Hashtable<String, ArrayList<MatchRow>> matches = new Hashtable<String, ArrayList<MatchRow>>();
	
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
	public void put(String key, MatchRow match) {
		ArrayList<MatchRow> existingMatches = matches.get(key);
		if (existingMatches == null) {
			existingMatches = new ArrayList<MatchRow>();
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
					ArrayList<MatchRow> matchArray = new ArrayList<MatchRow>();
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
	public ArrayList<MatchRow> getArrayList() {
		/* determining the size of the output */
		int size = 0;
		Enumeration<ArrayList<MatchRow>> e = matches.elements();
		while (e.hasMoreElements()) {
			size += e.nextElement().size();
		}
		
		/* create the output */
		ArrayList<MatchRow> out = new ArrayList<MatchRow>(size);
		e = matches.elements();
		while (e.hasMoreElements()) {
			out.addAll(e.nextElement());
		}
		
		return out;
	}
	
	public Hashtable<String, ArrayList<MatchRow>> getHashtable() {
		return matches;
	}
	
	public ArrayList<MatchRow> get(String key) {
		return matches.get(key);
	}

	
}
