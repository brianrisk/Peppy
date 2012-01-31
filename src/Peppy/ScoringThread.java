package Peppy;
import java.util.ArrayList;

import Math.MassError;
import Math.MathFunctions;


public class ScoringThread implements Runnable {
	
	/* this holds the full list of peptides (sorted by mass), though we will only be using a section for the given spectrum */
	ArrayList<Peptide> peptides;
	
	/* the spectrum for which we will be searching */
	Spectrum spectrum;
	
	/* our commanding server; results get reported back to this and new spectra to search come from here */
	ScoringThreadServer scoringThreadServer;
	
	/* due to the imperfections of my binary search, we need this extra margin */
	private static int extraMargin = 8;
	
	/**
	 * @param peptides
	 * @param spectrum
	 */
	public ScoringThread(Spectrum spectrum, ArrayList<Peptide> peptides, ScoringThreadServer scoringThreadServer) {
		this.spectrum = spectrum;
		this.peptides = peptides;
		this.scoringThreadServer = scoringThreadServer;
	}
	

	public void run() {
		
		while (spectrum != null) {
	
			//find the first index of the peptide with mass greater than lowestPeptideMassToConsider
			double lowestPeptideMassToConsider = spectrum.getMass() - MassError.getDaltonError(Properties.precursorTolerance, spectrum.getMass());
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
			double highestPeptideMassToConsider = spectrum.getMass() + MassError.getDaltonError(Properties.precursorTolerance, spectrum.getMass());
			if (Properties.searchModifications) {
				/* ditto above */
				highestPeptideMassToConsider -= Properties.modificationLowerBound;
			}
			int lastPeptideIndex = MathFunctions.findFirstIndexGreater(peptides, highestPeptideMassToConsider);
			lastPeptideIndex += extraMargin;
			if (lastPeptideIndex >= peptides.size()) lastPeptideIndex = peptides.size() - 1;
			
			/* set up where to hold the matches */
			int peptideCount = lastPeptideIndex - firstPeptideIndex;
			ArrayList<Match> matchesForEValueCalculation = new ArrayList<Match>();
			ArrayList<Match> topMatches = new ArrayList<Match>(peptideCount);
						
			//examine only peptides in our designated mass range
			for (int peptideIndex = firstPeptideIndex; peptideIndex < lastPeptideIndex; peptideIndex++) {
				Peptide peptide = peptides.get(peptideIndex);
				
				Match match = Properties.matchConstructor.createMatch(spectrum, peptide);
				
				/* add the match to our E value calculations list*/
				if (Properties.calculateEValues) {
					matchesForEValueCalculation.add(match);
				}
						
				/* only add the match if it is decent */
				if (match.getIMP() <= Properties.maxIMP) {
					topMatches.add(match);
				}
				
			}
			
			/* keep track of scores for e values */
			if (Properties.calculateEValues) {
				spectrum.getEValueCalculator().addScores(matchesForEValueCalculation);
			}

			/* return results, get new task */
			spectrum = scoringThreadServer.getNextSpectrum(topMatches);

		}
	}

	
	

}
