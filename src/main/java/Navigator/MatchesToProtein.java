package Navigator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

import Experimental.ProteinOverlap;


/**
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class MatchesToProtein extends MatchesTo  {
	
	
	private Hashtable<String, Match> bestPeptides = new Hashtable<String, Match>();
	
	private int uniquePeptideCount = 0;
	
	private double maxUniquePeptideScore = 0;
	
	private ArrayList<ProteinOverlap> proteinOverlaps;
	
	private boolean isEntirelyFoundInOtherProtein = false;
	
	private String chromosome = "";
	

	public MatchesToProtein(String name) {
		super(name);
	}

	public void addMatch(Match match) {
		super.addMatch(match);
		
		Match bestPeptideMatch = bestPeptides.get(match.getString("peptideSequence"));
		if (bestPeptideMatch == null) {
			bestPeptides.put(match.getString("peptideSequence"), match);
			if (match.getInt("RankCount") == 1) {
				uniquePeptideCount++;
				if (match.getScore() > maxUniquePeptideScore) maxUniquePeptideScore = match.getScore();
			}
		} else {
			if (bestPeptideMatch.getScore() < match.getScore()) {
				bestPeptides.put(match.getString("peptideSequence"), match);
			}
		}
		
	}
	

	public ArrayList<ProteinOverlap> getProteinOverlaps() {
		return proteinOverlaps;
	}

	public boolean hasUniqueMatch() {
		return (uniquePeptideCount > 0);
	}
	
	public int getUniquePeptideCount() {
		return uniquePeptideCount;
	}


	public double getMaxUniquePeptideScore() {
		return maxUniquePeptideScore;
	}

	public Hashtable<String, Match> getBestPeptides() {
		return bestPeptides;
	}

	
	/**
	 *  get proteins that have overlapping peptides as this protein 
	 */
	public void determineProteinOverlaps(Hashtable<String, ArrayList<MatchesToProtein>> proteinsForPeptide) {
		
		/* this is a hash so that we only have one entry for every protein that shares a peptide with this protein */
		Hashtable<String, MatchesToProtein> overLappingProteins = new Hashtable<String, MatchesToProtein>();
		
		/* now we get all the individual peptides in this protein */
		ArrayList<String> peptidesInThisProtein = new ArrayList<String>(getBestPeptides().keySet());
		
		/* we look up each peptide and find if it is in other proteins */
		for (String peptide: peptidesInThisProtein) {
			ArrayList<MatchesToProtein> proteinsWithThisPeptide = proteinsForPeptide.get(peptide);
			for (MatchesToProtein relatedProtein: proteinsWithThisPeptide) {
				
				/* skip if the protein is this protein */
				if (relatedProtein.getName().equals(getName())) continue;
				
				/* add the protein to our list of overlapping proteins */
				overLappingProteins.put(relatedProtein.getName(), relatedProtein);
			}
		}
		
		/* convert the overlapping proteins hash to an easier array list */
		ArrayList<MatchesToProtein> relatedProteins = new ArrayList<MatchesToProtein>(overLappingProteins.values());
		
		/* create some protein overlap objects and fill this list */
		proteinOverlaps = new ArrayList<ProteinOverlap>(relatedProteins.size());
		for (MatchesToProtein proteinB: relatedProteins) {
			ProteinOverlap overlap = new ProteinOverlap(this, proteinB);
			if (overlap.getOnlyACount() == 0) isEntirelyFoundInOtherProtein = true;
			proteinOverlaps.add(overlap);
		}
		Collections.sort(proteinOverlaps);
	}

	public boolean isEntirelyFoundInOtherProtein() {
		return isEntirelyFoundInOtherProtein;
	}



}
