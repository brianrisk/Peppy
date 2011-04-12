package Peppy;
import java.util.ArrayList;
import java.util.Collections;

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
	
			//find the first index of the peptide with mass greater than lowestPeptideMassToConsider
			int firstPeptideIndex;
			if (Properties.scoringMethodName.equals("Peppy.Match_IMP_MultiMod")) {
				firstPeptideIndex = 0;
			} else {
				double lowestPeptideMassToConsider = spectrum.getMass() - Properties.spectrumToPeptideMassError;
				firstPeptideIndex = MathFunctions.findFirstIndexGreater(peptides, lowestPeptideMassToConsider);
				firstPeptideIndex -= 8;
				if (firstPeptideIndex < 0) firstPeptideIndex = 0;
			}
			
			
			//find the last index, compensate for rounding error
			double highestPeptideMassToConsider = spectrum.getMass() + Properties.spectrumToPeptideMassError;
			int lastPeptideIndex = MathFunctions.findFirstIndexGreater(peptides, highestPeptideMassToConsider);
			lastPeptideIndex += 8;
			if (lastPeptideIndex >= peptides.size()) lastPeptideIndex = peptides.size() - 1;
			
			
			//examine only peptides in our designated mass range
			for (int peptideIndex = firstPeptideIndex; peptideIndex < lastPeptideIndex; peptideIndex++) {
				Peptide peptide = peptides.get(peptideIndex);
				
				Match match = Properties.matchConstructor.createMatch(spectrum, peptide);

				if (match.getScore() == 0.0) continue;
				
				matchesForOneSpectrum.add(match);
			}
			
			//collect the top maximumNumberOfMatchesForASpectrum
			Match.setSortParameter(Match.SORT_BY_SCORE);
			Collections.sort(matchesForOneSpectrum);
			ArrayList<Match> topMatches = new ArrayList<Match>();
			int max = Properties.maximumNumberOfMatchesForASpectrum;
			if (matchesForOneSpectrum.size() < max) max = matchesForOneSpectrum.size();
			for (int i = 0; i < max; i++) {
				topMatches.add(matchesForOneSpectrum.get(i));
			}
			
			//We may want to reduce the number of duplicate matches
			if (Properties.reduceDuplicateMatches && max > 1) {
				for (int i = 0; i < topMatches.size(); i++) {
					Match matchA = topMatches.get(i);
					for (int j = i + 1; j < topMatches.size(); j++) {
						Match matchB = topMatches.get(j);
						if (matchA.equals(matchB)) {
							topMatches.remove(j);
							j--;
						}
					}
				}
			}

			
			//assign E values to top Matches:
			if (matchesForOneSpectrum.size() != 0) {
				spectrum.calculateEValues(matchesForOneSpectrum, topMatches);
			}
			
			//return results, get new task
			spectrum = scoringThreadServer.getNextSpectrum(topMatches);
		}
	}

	
	

}
