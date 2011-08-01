package Peppy;
import java.util.ArrayList;

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
			if (Properties.scoringMethodName.equals("Peppy.Match_IMP_MultiMod")) {
				firstPeptideIndex = 0;
			} else {
				double lowestPeptideMassToConsider = spectrum.getMass() - Properties.precursorTolerance;
				firstPeptideIndex = MathFunctions.findFirstIndexGreater(peptides, lowestPeptideMassToConsider);
				firstPeptideIndex -= 8;
				if (firstPeptideIndex < 0) firstPeptideIndex = 0;
			}
			
			
			//find the last index, compensate for rounding error
			double highestPeptideMassToConsider = spectrum.getMass() + Properties.precursorTolerance;
			int lastPeptideIndex = MathFunctions.findFirstIndexGreater(peptides, highestPeptideMassToConsider);
			lastPeptideIndex += 8;
			if (lastPeptideIndex >= peptides.size()) lastPeptideIndex = peptides.size() - 1;
			
			//examine only peptides in our designated mass range
			for (int peptideIndex = firstPeptideIndex; peptideIndex < lastPeptideIndex; peptideIndex++) {
				Peptide peptide = peptides.get(peptideIndex);
				
				Match match = Properties.matchConstructor.createMatch(spectrum, peptide);
				
				/* add the match we find */
				matchesForOneSpectrum.add(match);
						
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
