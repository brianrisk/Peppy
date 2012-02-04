package Navigator;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;

import Peppy.Match;


public class GroupingPeptide {
	
	Hashtable<String, ArrayList<Match>> peptides = new Hashtable<String, ArrayList<Match>>();
	
	/**
	 * Constructor
	 */
	public GroupingPeptide() {
		
	}
	
	/**
	 * Adding a match to the group
	 * @param match
	 */
	public void addMatch(Match match) {
		ArrayList<Match> group = peptides.get(match.getPeptide().getAcidSequenceString());
		
		/* if there is no entry for that peptide, make one */
		if (group == null) {
			group = new ArrayList<Match>();
			group.add(match);
			peptides.put(match.getPeptide().getAcidSequenceString(), group);
		} 
		
		/* else, retrieve list and add to it */
		else {
			group.add(match);
		}
	}
	
	
	/**
	 * Sorting the matches for each peptide
	 */
	public void sort() {
		Enumeration<ArrayList<Match>> values = peptides.elements();
		while (values.hasMoreElements()) {
			Collections.sort(values.nextElement());
		}
	}

}
