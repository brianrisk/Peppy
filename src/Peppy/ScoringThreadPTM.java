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
			
			ArrayList<Match_IMP_VariMod> matchesForOneSpectrum = new ArrayList<Match_IMP_VariMod>();
	
			//find the first index of the peptide with mass greater than lowestPeptideMassToConsider
			double lowestSpectrumMassToConsider = peptide.getMass() - 40;
			int firstSpectrumIndex = MathFunctions.findFirstIndexGreater(spectra, lowestSpectrumMassToConsider);
			firstSpectrumIndex -= 8;
			if (firstSpectrumIndex < 0) firstSpectrumIndex = 0;
			
			//find the last index, compensate for rounding error
			double highestSpectrumMassToConsider = peptide.getMass() + 550;
			int lastSpectrumIndex = MathFunctions.findFirstIndexGreater(spectra, highestSpectrumMassToConsider);
			lastSpectrumIndex += 8;
			if (lastSpectrumIndex >= spectra.size()) lastSpectrumIndex = spectra.size() - 1;
			
			//examine only peptides in our designated mass range
			for (int spectrumIndex = firstSpectrumIndex; spectrumIndex < lastSpectrumIndex; spectrumIndex++) {
				Spectrum spectrum = spectra.get(spectrumIndex);
				Match_IMP_VariMod match = new Match_IMP_VariMod(spectrum, peptide);
				if (match.getScore() == 1.0) {
					continue;
				}
				matchesForOneSpectrum.add(match);
			}
			
			//collect the top maximumNumberOfMatchesForASpectrum
			Match_IMP_VariMod.setSortParameter(Match_IMP_VariMod.SORT_BY_SCORE);
			Collections.sort(matchesForOneSpectrum);
			ArrayList<Match_IMP_VariMod> topMatches = new ArrayList<Match_IMP_VariMod>();
			int max = Properties.maximumNumberOfMatchesForASpectrum;
			if (matchesForOneSpectrum.size() < max) max = matchesForOneSpectrum.size();
			for (int i = 0; i < max; i++) {
				topMatches.add(matchesForOneSpectrum.get(i));
			}
			
			//We may want to reduce the number of duplicate matches
			if (Properties.reduceDuplicateMatches && max > 1) {
				for (int i = 0; i < topMatches.size(); i++) {
					Match_IMP_VariMod matchA = topMatches.get(i);
					for (int j = i + 1; j < topMatches.size(); j++) {
						Match_IMP_VariMod matchB = topMatches.get(j);
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
