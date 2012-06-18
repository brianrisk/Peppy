package Peppy;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;

import Math.MathFunctions;
import Reports.HTMLPageRegions;

/**
 * Copyright 2012, Brian Risk
 * Released under the Netscape Public License
 * 
 * @author Brian Risk
 *
 */
public class Regions {
	
	private int maxLength = 1024;
	ArrayList<Match> matches;
	ArrayList<Sequence> sequences;
	
	/* where we store the interesting regions */
	ArrayList<Region> regions;
	
	
	public Regions(ArrayList<Match> matches, ArrayList<Sequence> sequences) {
		this.matches = matches;
		this.sequences = sequences;
		
		/* Find interesting regions */
		if (Properties.isSequenceFileDNA) {
			
			/* set sort parameter to sort by locus */
			Match.setSortParameter(Match.SORT_BY_LOCUS);
			
			/* finds forwards strand and reverse strand regions */
			ArrayList<Region> forwardRegions = findRegions(true);
			ArrayList<Region> reverseRegions = findRegions(false);
			regions = new ArrayList<Region>(forwardRegions.size() + reverseRegions.size());
			regions.addAll(forwardRegions);
			regions.addAll(reverseRegions);
			Collections.sort(regions);
			
		}
	}
	
	/**
	 * finds regions on forwards or reverse strand
	 * @param isForward
	 */
	private ArrayList<Region> findRegions(boolean isForward) {
		
		ArrayList<Region> returnedRegions = new ArrayList<Region>();
		
		/* get all interesting regions from each sequence */
		for (Sequence sequence: sequences) {
			ArrayList<Match> sequenceMatches = Matches.getMatchesWithSequence((Sequence_DNA) sequence, matches);
			Collections.sort(sequenceMatches);
			
			/* tracks which regions we are still building 
			 * initializes with every new sequence so that
			 * we don't add matches from one sequence to 
			 * regions in another */
			int activeRegionIndex = returnedRegions.size();

			int previousLocus = -1;
			int locus = -1;
			for (Match match: sequenceMatches) {
				if (match.getPeptide().isForward() != isForward) continue;
				locus = match.getPeptide().getStartIndex();
				
				/* if we've got a different locus, that means some regions may be finished */
				if (locus != previousLocus) {
					
					/* create a new region */
					returnedRegions.add(new Region(locus, maxLength, sequence, isForward));
					
					/* find our new index where we are still building */
					for (int i = activeRegionIndex; i < returnedRegions.size(); i++) {
						Region region = returnedRegions.get(i);
						if (region.getStopLocation() > locus) {
							activeRegionIndex = i;
							break;
						}
					}
					
					/* update previous locus */
					previousLocus = locus;
				}
				
				/* add match to all active regions */
				for (int i = activeRegionIndex; i < returnedRegions.size(); i++) {
					Region region = returnedRegions.get(i);
					region.addMatch(match);
				}
			}
		}
			
		
		/* find best non-overlapping regions */
		Collections.sort(returnedRegions);
		if (returnedRegions.size() > 1) {
			Region regionA, regionB;
			for (int i = 0; i < returnedRegions.size() - 1; i++) {
				regionA = returnedRegions.get(i);
				if (regionA.isFlagged()) continue;
				for (int j = i + 1; j < returnedRegions.size(); j++) {
					regionB = returnedRegions.get(j);
					if (regionA.isOverlapping(regionB)) {
						regionB.flag();
					}
				}
			}
		}
		
		/* pull out the non-flagged regions */
		ArrayList<Region> nonOverlapping = new ArrayList<Region>();
		for (Region region: returnedRegions) {
			if (region.isUnFlagged()) {
				nonOverlapping.add(region);
			}
		}
		returnedRegions = nonOverlapping;
		
		/* calculate scores */
		for (Region region: returnedRegions) {
			region.calculateScore();
		}
		
		/* calculate P values */
		double pValue;
		
		/* this will actually compensate for the fact that we will have less
		 * regions in reverse searches due to E value cutoffs.  If we assumed the
		 * same amount of regions (or maxlength / genome size) then it would appear
		 * that probabilities for the reverse search were lower than they actually are
		 */
		double probability = 1.0 / (double) returnedRegions.size();
//		double probability = (double) maxLength / 3000000000.0; // region divided by length of genome 
		for (Region region: returnedRegions) {
//			eValue = MathFunctions.approximateBinomialProbability(matches.size(), region.getNumberOfMatches(), probability);
			pValue = MathFunctions.approximateNegativeLog10OfBinomialProbability(matches.size(), region.getNumberOfMatches(), probability);
			region.setPValue(pValue);
		}
		
		/* I know the regions should already be sorted by now, but this just makes me feel better  */
		Collections.sort(returnedRegions);
		
		return returnedRegions;
	}
	
	
	
	public void createReport(File reportDir) {
			
			/* set up the regions folders */
			File regionsHTMLFolder = new File(reportDir, "regions html");
			regionsHTMLFolder.mkdirs();
			
			/* determine how many regions we are reporting */
			int max = Properties.numberOfRegionsToKeep;
			if (regions.size() < max) max = regions.size();
			
			/* trim down regions to this amount */
			ArrayList<Region> lessRegions = new ArrayList<Region>(max);
			for (int i = 0; i < max; i++) {
				lessRegions.add(regions.get(i));
			}
			regions = lessRegions;
			
			/* set up the main regions HTML file */
			HTMLPageRegions indexPage = new HTMLPageRegions(new File(regionsHTMLFolder, "index.html"), regions);
			indexPage.makePage();

	}
	
	
	
	public ArrayList<Region> getRegions() {
		return regions;
	}

	public void clearRegions() {
		for (Region region: regions) {
			region.clearRegion();
		}
		regions.clear();
	}

}
