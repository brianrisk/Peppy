package Peppy;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import Math.MassError;
import Math.MathFunctions;

/**
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class ScoringThread implements Runnable {
	
	/* this holds the full list of peptides (sorted by mass), though we will only be using a section for the given spectrum */
	ArrayList<Peptide> peptides;
	
	/* Where we will be holding our matches */
	MatchesSpectrum matchesSpectrum;
	
	/* our commanding server; results get reported back to this and new spectra to search come from here */
	ScoringServer scoringServer;
	
	/* due to the imperfections of my binary search, we need this extra margin */
	private static int extraMargin = 0;
	
	/**
	 * @param peptides
	 * @param spectrum
	 */
	public ScoringThread(MatchesSpectrum matchesSpectrum, ArrayList<Peptide> peptides, ScoringServer scoringServer) {
		this.matchesSpectrum = matchesSpectrum;
		this.peptides = peptides;
		this.scoringServer = scoringServer;
	}
	
	
//	public static void main(String args[] ) {
//		Random random = new Random();
//		ArrayList<Double> numbers = new ArrayList<Double>();
//		for (int i = 0; i < 100; i++) {
//			numbers.add(random.nextDouble());
//		}
//		numbers.add(0.5);
//		Collections.sort(numbers);
//		int location = Collections.binarySearch(numbers, 0.5);
//		if (location < 0) {
//			location++;
//			location = Math.abs(location);
//		}
//		U.p(numbers.get(location - 1));
//		U.p(numbers.get(location));
//		U.p(numbers.get(location + 1));
//	}
	
	public void run() {
		
		while (matchesSpectrum != null) {
			if (Properties.useHydrogenOffset) {
				for (int index = 0; index < 3; index++) {
					searchPeptideRange(index);
				}
			} else {
				searchPeptideRange(0);
			}
	
			/* return results, get new spectrum to process.  
			 * If no more spectra, gets null and the loop exits */
			matchesSpectrum = scoringServer.getNextSpectrumMatches();
	
		}
	}


	private void searchPeptideRange(int hydrogenOffset) {
		
		//find the first index of the peptide with mass greater than lowestPeptideMassToConsider
		double lowestPeptideMassToConsider = matchesSpectrum.getSpectrum().getMass() - MassError.getDaltonError(Properties.precursorTolerance, matchesSpectrum.getSpectrum().getMass());
		lowestPeptideMassToConsider -= Definitions.HYDROGEN_MONO * hydrogenOffset;
		if (Properties.searchModifications) {
			/* I know subtracting the upper bound seems backwards, but since a 
			 * modification on the peptide makes the spectrum heavier, this is the
			 * order things should be*/
			lowestPeptideMassToConsider -= Properties.modificationUpperBound;
		}
		int firstPeptideIndex = MathFunctions.findFirstIndexGreater(peptides, lowestPeptideMassToConsider);
		firstPeptideIndex -= extraMargin;
		if (firstPeptideIndex < 0) firstPeptideIndex = 0;
		
		
		//find the last index, compensate for rounding error
		double highestPeptideMassToConsider = matchesSpectrum.getSpectrum().getMass() + MassError.getDaltonError(Properties.precursorTolerance, matchesSpectrum.getSpectrum().getMass());
		highestPeptideMassToConsider -= Definitions.HYDROGEN_MONO * hydrogenOffset;
		if (Properties.searchModifications) {
			/* ditto above */
			highestPeptideMassToConsider -= Properties.modificationLowerBound;
		}
		int lastPeptideIndex = MathFunctions.findFirstIndexGreater(peptides, highestPeptideMassToConsider);
		lastPeptideIndex += extraMargin;
		if (lastPeptideIndex >= peptides.size()) lastPeptideIndex = peptides.size() - 1;
					
		/* examine only peptides in our designated mass range */
		for (int peptideIndex = firstPeptideIndex; peptideIndex < lastPeptideIndex; peptideIndex++) {
			
			Match match = Properties.matchConstructor.createMatch(matchesSpectrum, peptides.get(peptideIndex));
			matchesSpectrum.addMatch(match);
			
		}
	}
	

	private int findIndex(double mass) {
		int index = Collections.binarySearch(peptides, new Peptide(mass));
		if (index < 0) index = - index - 1;
		if (index < 0) index = 0;
		if (index >= peptides.size()) index = peptides.size() - 1;
		return index;
	}

	
	

}
