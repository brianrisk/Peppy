package Peppy;
import java.util.ArrayList;

import Math.MassError;
import Math.MathFunctions;


public class ScoringThread implements Runnable {
	
	ArrayList<Peptide> peptides;
	Spectrum spectrum;
	ScoringThreadServer scoringThreadServer;
	
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
			
			ArrayList<Match> matchesForOneSpectrum = new ArrayList<Match>();
			ArrayList<Match> topMatches = new ArrayList<Match>();
	
			//find the first index of the peptide with mass greater than lowestPeptideMassToConsider
			int firstPeptideIndex;
			double lowestPeptideMassToConsider = spectrum.getMass() - MassError.getDaltonError(Properties.precursorTolerance, spectrum.getMass());
			if (Properties.searchModifications) {
				/* I know subtracting the upper bound seems backwards, but since a 
				 * modification on the peptide makes the spectrum heavier, this is the
				 * order things should be*/
				lowestPeptideMassToConsider -= Properties.modificationUpperBound;
			}
			firstPeptideIndex = MathFunctions.findFirstIndexGreater(peptides, lowestPeptideMassToConsider);
			firstPeptideIndex -= 8;
			if (firstPeptideIndex < 0) firstPeptideIndex = 0;
			
			
			//find the last index, compensate for rounding error
			double highestPeptideMassToConsider = spectrum.getMass() + MassError.getDaltonError(Properties.precursorTolerance, spectrum.getMass());
			if (Properties.searchModifications) {
				/* ditto above */
				highestPeptideMassToConsider -= Properties.modificationLowerBound;
			}
			int lastPeptideIndex = MathFunctions.findFirstIndexGreater(peptides, highestPeptideMassToConsider);
			lastPeptideIndex += 8;
			if (lastPeptideIndex >= peptides.size()) lastPeptideIndex = peptides.size() - 1;
			
			//examine only peptides in our designated mass range
			for (int peptideIndex = firstPeptideIndex; peptideIndex < lastPeptideIndex; peptideIndex++) {
				Peptide peptide = peptides.get(peptideIndex);
				
				Match match = Properties.matchConstructor.createMatch(spectrum, peptide);
				
				/* add the match we find */
				if (match.getIMP() <= 1.0E-8) {
					matchesForOneSpectrum.add(match);
				}
						
				/* only add the match if it is decent */
				if (match.getIMP() <= Properties.maxIMP) {
					topMatches.add(match);
				}
				
			}
			
			/* keep track of scores for e values */
			spectrum.getEValueCalculator().addScores(matchesForOneSpectrum);

			/* return results, get new task */
			spectrum = scoringThreadServer.getNextSpectrum(topMatches);

		}
	}

	
	

}
