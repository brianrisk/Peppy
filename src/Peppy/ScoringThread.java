package Peppy;
import java.util.ArrayList;
import java.util.Collections;


public class ScoringThread implements Runnable {
	
	ArrayList<Peptide> peptides;
	Spectrum spectrum;
	ScoringThreadServer scoringThreadServer;
	Sequence sequence;
	
	/**
	 * @param peptides
	 * @param spectrum
	 */
	public ScoringThread(Spectrum spectrum, ArrayList<Peptide> peptides, ScoringThreadServer scoringThreadServer, Sequence sequence) {
		this.spectrum = spectrum;
		this.peptides = peptides;
		this.scoringThreadServer = scoringThreadServer;
		this.sequence = sequence;
	}

	public void run() {
		
		while (spectrum != null) {
			
			ArrayList<Match> matchesForOneSpectrum = new ArrayList<Match>();
	
			//find the first index of the peptide with mass greater than lowestPeptideMassToConsider
			double lowestPeptideMassToConsider = spectrum.getPrecursorMass() - Properties.spectrumToPeptideMassError;
			int firstPeptideIndex = findFirstIndexWithGreaterMass(peptides, lowestPeptideMassToConsider);
			
			//find the first index of the peptide with mass greater than highestPeptideMassToConsider
			double highestPeptideMassToConsider = spectrum.getPrecursorMass() + Properties.spectrumToPeptideMassError;
			int lastPeptideIndex = findFirstIndexWithGreaterMass(peptides, highestPeptideMassToConsider);
			
			//examine only peptides in our designated mass range
			for (int peptideIndex = firstPeptideIndex; peptideIndex < lastPeptideIndex; peptideIndex++) {
				Peptide peptide = peptides.get(peptideIndex);
				Match match = new Match(spectrum, peptide, sequence);
				if (match.getScore() == 0.0) {
					continue;
				}
				matchesForOneSpectrum.add(match);
			}
			
//			if (spectrum.getFile().getName().equals("T10707_Well_H13_1768.77_19185.mgf..pkl")) {
//				U.p("first index: " + firstPeptideIndex);
//				U.p("last index: " + lastPeptideIndex);
//				U.p("precursor = " + spectrum.getPrecursorMass());
//				U.p ("lowestPeptideMassToConsider: " + lowestPeptideMassToConsider);
//				U.p ("highestPeptideMassToConsider: " + highestPeptideMassToConsider);
//				U.p("peptide mass = " + (new Peptide("VVSMDENFHPLNELIPLVYIQDPK")).getMass());
//			}
			
			//collect the top maximumNumberOfMatchesForASpectrum
			Match.setSortParameter(Match.SORT_BY_DEFAULT);
			Collections.sort(matchesForOneSpectrum);
			ArrayList<Match> topMatches = new ArrayList<Match>();
			int max = Properties.maximumNumberOfMatchesForASpectrum;
			if (matchesForOneSpectrum.size() < max) max = matchesForOneSpectrum.size();
			
				
//				//make sure adjacent matches aren't repeated
//				Match topMatch = matchesForOneSpectrum.get(0);
//				Match presentMatch;
//				Match previousMatch = matchesForOneSpectrum.get(0);
//				topMatches.add(matchesForOneSpectrum.get(0));
//				for (int i = 1; i < max; i++) {
//					presentMatch = matchesForOneSpectrum.get(i);
//					if (!presentMatch.equals(previousMatch)) {
//						topMatches.add(presentMatch);
//					}else {
//						//since we're skipping one we are increasing the upper bound	
//						if (max < matchesForOneSpectrum.size()) max++;
//					}
//					previousMatch = presentMatch;
//				}

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

			//set rank -- NOTE: this might not be true rank if multiple chromosomes, multiple digestion windows
			for (int i = 0; i < topMatches.size(); i++) {
				Match topMatch = topMatches.get(i);
				topMatch.setRank(i);
			}
			
			
			//assign E values to top Matches:
			if (matchesForOneSpectrum.size() != 0) {
				spectrum.calculateEValues(matchesForOneSpectrum, topMatches);
			}
			
			//return results, get new task
			spectrum = scoringThreadServer.getNextSpectrum(topMatches);
		}
	}

	
	/**
	 * Boolean search to locate the first peptide in the SORTED list of peptides that has
	 * a mass greater than the "mass" parameter.
	 * @param peptides
	 * @param mass
	 * @return
	 */
	public static int findFirstIndexWithGreaterMass(ArrayList<Peptide> peptides, double mass) {
		Peptide peptide;
		int index = peptides.size() / 2;
		int increment = index / 2;
		while (increment > 0) {
			peptide = peptides.get(index);
			if (peptide.getMass() > mass) {index -= increment;}
			else {index += increment;}
			increment /= 2;
		}
		return index;
	}

}
