package Peppy;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;

import Math.EValueCalculator;
import Math.MathFunctions;
import Reports.CommonMatchSearches;
import Reports.HTMLPageRegions;
import Reports.HistogramVisualizer;
import Utilities.U;

public class Regions {
	
	/* where we store the interesting regions */
	ArrayList<Region> regions = new ArrayList<Region>();
	
	/* E values */
	EValueCalculator evc = new EValueCalculator();
	
	
	public Regions(ArrayList<Match> matches, ArrayList<Sequence> sequences, ArrayList<Spectrum> spectra) {
		/* Find interesting regions */
		if (Properties.isSequenceFileDNA) {
			
			/* set sort parameter to sort by locus */
			Match.setSortParameter(Match.SORT_BY_LOCUS);
			
			
			/* get all interesting regions from each sequence */
			for (Sequence sequence: sequences) {
				ArrayList<Match> sequenceMatches = CommonMatchSearches.getMatchesWithSequence((Sequence_DNA) sequence, matches);
				Collections.sort(sequenceMatches);
				
				/* tracks which regions we are still building 
				 * initializes with every new sequence so that
				 * we don't add matches from one sequence to 
				 * regions in another */
				int activeRegionIndex = regions.size();

				int previousLocus = -1;
				int locus = -1;
				for (Match match: sequenceMatches) {
					locus = match.getPeptide().getStartIndex();
					
					/* if we've got a different locus, that means some regions may be finished */
					if (locus != previousLocus) {
						
						/* create a new region */
						regions.add(new Region(locus, 1024, sequence));
						
						/* find our new index where we are still building */
						for (int i = activeRegionIndex; i < regions.size(); i++) {
							Region region = regions.get(i);
							if (region.getStopLocation() > locus) {
								activeRegionIndex = i;
								break;
							}
						}
						
						/* update previous locus */
						previousLocus = locus;
					}
					
					/* add match to all active regions */
					for (int i = activeRegionIndex; i < regions.size(); i++) {
						Region region = regions.get(i);
						region.addMatch(match);
					}
				}
			}
				
			
			/* find best non-overlapping regions */
			Collections.sort(regions);
			if (regions.size() > 1) {
				Region regionA, regionB;
				for (int i = 0; i < regions.size() - 1; i++) {
					regionA = regions.get(i);
					if (regionA.isFlagged()) continue;
					for (int j = i + 1; j < regions.size(); j++) {
						regionB = regions.get(j);
						if (regionA.isOverlapping(regionB)) {
							regionB.flag();
						}
					}
				}
			}
			
			/* pull out the non-flagged regions */
			ArrayList<Region> nonOverlapping = new ArrayList<Region>();
			for (Region region: regions) {
				if (region.isUnFlagged()) {
					nonOverlapping.add(region);
				}
			}
			regions = nonOverlapping;
			
			/* calculate scores */
			for (Region region: regions) {
				region.calculateScore();
			}
			
			/* calculate E values */
//			evc.addScores(regions);
//			evc.calculateHistogramProperties();
//			for (Region region: regions) {
//				region.setEValue(evc.calculateEValueOfScore(region.getScore()));
//			}
			double eValue;
			double probability = 1.0 / (double) regions.size();
			for (Region region: regions) {
				eValue = MathFunctions.approximateBinomialProbability(spectra.size(), region.getNumberOfMatches(), probability);
				region.setEValue(eValue);
			}

			
		}
	}
	
	public void createReport(File reportDir) {
		try {
			
			
			
			/* set up the regions folders */
			File regionsTextFolder = new File(reportDir, "regions");
			regionsTextFolder.mkdirs();
			File regionsHTMLFolder = new File(reportDir, "regions html");
			regionsHTMLFolder.mkdirs();
			
			/* determine how many regions we are reporting */
			int max = 1000;
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
			
			
			for (int i = 0; i < max; i ++ ) {
				Region region = regions.get(i);
				String spacer = "";
				if (i < 1000) {
					spacer += "0";
					if (i < 100) {
						spacer += "0";
						if (i < 10) {
							spacer += "0";
						}
					}
				}
				int neatE = (int) Math.round(-Math.log10(region.getEValue()));
				PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(regionsTextFolder, "region " + spacer + i + " " + neatE + ".txt"))));
				ArrayList<Match> regionMatches = region.getMatches();
				for (Match match: regionMatches) {
					pw.println(match);
				}
				pw.flush();
				pw.close();
			}
			
			/* make the histogram */
			File histogramFile = new File(regionsHTMLFolder, "histogram.jpg");
			HistogramVisualizer.drawHistogram(evc.getSmoothedHistogram(), 300, 300, histogramFile);
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void clearRegions() {
		for (Region region: regions) {
			region.clearRegion();
		}
		regions.clear();
	}

}
