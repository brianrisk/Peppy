package Peppy;
import java.util.ArrayList;
import java.util.Collections;

import Math.MathFunctions;


public class ScoringThreadPTM implements Runnable {
	
	ArrayList<Spectrum> spectra;
	Peptide peptide;
	ScoringThreadServerPTM scoringThreadServerPTM;
	
	/**
	 * @param peptides
	 * @param spectrum
	 */
	public ScoringThreadPTM(Peptide peptide, ArrayList<Spectrum> spectra, ScoringThreadServerPTM scoringThreadServerPTM) {
		this.spectra = spectra;
		this.peptide = peptide;
		this.scoringThreadServerPTM = scoringThreadServerPTM;
	}

	/*
	 * (non-Javadoc)
	 * CHANGES TO NOTE:
	 * match.getScore() == 1.0
	 * Match now is MatchPTM
	 * @see java.lang.Runnable#run()
	 */
	public void run() {
		
		while (peptide != null) {
			
			ArrayList<MatchPTM> matchesForOneSpectrum = new ArrayList<MatchPTM>();
	
			//find the first index of the peptide with mass greater than lowestPeptideMassToConsider
			double lowestSpectrumMassToConsider = peptide.getMass() - 20;
			int firstSpectrumIndex = MathFunctions.findFirstIndexGreater(spectra, lowestSpectrumMassToConsider);
			firstSpectrumIndex -= 8;
			if (firstSpectrumIndex < 0) firstSpectrumIndex = 0;
			
			//find the last index, compensate for rounding error
			double highestSpectrumMassToConsider = peptide.getMass() + 200;
			int lastSpectrumIndex = MathFunctions.findFirstIndexGreater(spectra, highestSpectrumMassToConsider);
			lastSpectrumIndex += 8;
			if (lastSpectrumIndex >= spectra.size()) lastSpectrumIndex = spectra.size() - 1;
			
			//examine only peptides in our designated mass range
			for (int spectrumIndex = firstSpectrumIndex; spectrumIndex < lastSpectrumIndex; spectrumIndex++) {
				Spectrum spectrum = spectra.get(spectrumIndex);
				MatchPTM match = new MatchPTM(spectrum, peptide);
				if (match.getScore() == 1.0) {
					continue;
				}
				matchesForOneSpectrum.add(match);
			}
			
			//collect the top maximumNumberOfMatchesForASpectrum
			MatchPTM.setSortParameter(MatchPTM.SORT_BY_SCORE);
			Collections.sort(matchesForOneSpectrum);
			ArrayList<MatchPTM> topMatches = new ArrayList<MatchPTM>();
			int max = Properties.maximumNumberOfMatchesForASpectrum;
			if (matchesForOneSpectrum.size() < max) max = matchesForOneSpectrum.size();
			for (int i = 0; i < max; i++) {
				topMatches.add(matchesForOneSpectrum.get(i));
			}
			
			//We may want to reduce the number of duplicate matches
			if (Properties.reduceDuplicateMatches && max > 1) {
				for (int i = 0; i < topMatches.size(); i++) {
					MatchPTM matchA = topMatches.get(i);
					for (int j = i + 1; j < topMatches.size(); j++) {
						MatchPTM matchB = topMatches.get(j);
						if (matchA.equals(matchB)) {
							topMatches.remove(j);
							j--;
						}
					}
				}
			}

			
			//assign E values to top Matches:
//			if (matchesForOneSpectrum.size() != 0) {
//				peptide.calculateEValues(matchesForOneSpectrum, topMatches);
//			}
			
			//return results, get new task
			peptide = scoringThreadServerPTM.getNextPeptide(topMatches);
		}
	}

	
	

}
