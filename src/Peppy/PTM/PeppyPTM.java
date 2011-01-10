package Peppy.PTM;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;

import Peppy.Match;
import Peppy.Peppy;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.ProteinDigestion;
import Peppy.RNA_Digestor;
import Peppy.RNA_Sequence;
import Peppy.ScoringThreadServer;
import Peppy.Sequence;
import Peppy.Spectrum;
import Reports.HTMLReporter;
import Reports.TextReporter;
import Utilities.U;


/**
 * A version of Peppy which is given a list of probable peptides and searches
 * for those peptides amongst the spectra assuming PTMs have taken place
 * @author Brian Risk
 *
 */
public class PeppyPTM {
	
	public static void main(String [] args) {
		Peppy.init(args);
		new PeppyPTM(args);
		U.p("done");
	}
	
	public PeppyPTM(String [] args) {
		U.startStopwatch();
		
		//Load our spectra
		U.p("loading spectra...");
		ArrayList<Spectrum> spectra = Spectrum.loadSpectra();
		U.p("loaded " +spectra.size() + " spectra.");
		
		//load our list of peptides
		ArrayList<Peptide> peptides = new ArrayList<Peptide>();
		
		//initialize our ArrayList of matches
		ArrayList<Match> matches = new ArrayList<Match>();
		
		
		matches.addAll(getMatches(peptides, spectra));
		
		//create new report directory
		File reportDir = new File(Properties.reportDirectory, Properties.reportDirectoryTitle + " " + System.currentTimeMillis());
		
		U.p("creating text reports");
		TextReporter textReport = new TextReporter(matches, spectra, sequences, reportDir);
		textReport.generateFullReport();
		textReport.generatePropertiesFile();
		
		U.p();
		U.stopStopwatch();
	}
	

	
	public static ArrayList<Match> getMatches(Sequence sequence, ArrayList<Spectrum> spectra) {
			U.p("working on sequence: " +sequence.getSequenceFile().getName());
			ArrayList<Match> matches = new ArrayList<Match>() ;
			ArrayList<Peptide> peptides;
			if (Properties.isSequenceFileDNA) {
				peptides = sequence.extractMorePeptides();
				//continually extract peptides from the sequence until there aren't anymore
				while (peptides != null) {
					//This is where the bulk of the processing in long jobs takes
					ArrayList<Match> newMatches = (new ScoringThreadServer(peptides, spectra, sequence)).getMatches();
					//Possible to add only matches with a decent e value
					if (Properties.useEValueCutOff) {
						for (Match match: newMatches) {
							if (match.getEValue() <= Properties.eValueCutOff) matches.add(match);
						}
					} else {
						matches.addAll(newMatches);
					}
					//free up the memory of the old peptide arraylist
					peptides.clear();
					System.gc();
					peptides = sequence.extractMorePeptides();
				}
				sequence.clearNucleotideData();
				U.p("removing duplicate matches");
				removeDuplicateMatches(matches);
			} else {
				peptides = ProteinDigestion.getPeptidesFromDatabase(sequence.getSequenceFile());
				//This is where the bulk of the processing in long jobs takes
				ArrayList<Match> newMatches = (new ScoringThreadServer(peptides, spectra, sequence)).getMatches();
				//Add only matches with a decent e value
				if (Properties.useEValueCutOff) {
					for (Match match: newMatches) {
						if (match.getEValue() <= Properties.eValueCutOff) matches.add(match);
					}
				} else {
					matches.addAll(newMatches);
				}
			}
			U.p("assigning final match ranks");
			assignRankToMatches(matches);
			assignRepeatedPeptideCount(matches);

			return matches;
	}
	
	/**
	 * Gets matches where a list of peptides is already derived
	 * @param peptides
	 * @param spectra
	 * @param sequence
	 * @return
	 */
	public static ArrayList<Match> getMatches(ArrayList<Peptide> peptides, ArrayList<Spectrum> spectra, Sequence sequence) {
		ArrayList<Match> matches = new ArrayList<Match>() ;
		//This is where the bulk of the processing in long jobs takes
		ArrayList<Match> newMatches = (new ScoringThreadServer(peptides, spectra, sequence)).getMatches();
		//Add only matches with a decent e value
		if (Properties.useEValueCutOff) {
			for (Match match: newMatches) {
				if (match.getEValue() <= Properties.eValueCutOff) matches.add(match);
//				if (match.getScore() >= 40.0) matches.add(match);
			}
		} else {
			matches.addAll(newMatches);
		}
		
		U.p("removing duplicate matches");
		removeDuplicateMatches(matches);
		
		U.p("assigning final match ranks");
		assignRankToMatches(matches);
		
		return matches;
	}
	
		
	public static void assignRankToMatches(ArrayList<Match> matches) {
		//first error check
		if (matches.size() == 0) return;
		Match.setSortParameter(Match.SORT_BY_SPECTRUM_ID_THEN_SCORE);
		Collections.sort(matches);
		Match match = matches.get(0);
		Match previousMatch = matches.get(0);
		//set for the first
		match.setRank(1);
		int rank = 1;
		for (int i = 1; i < matches.size(); i++) {
			//see if these are matches for a different spectrum
			match = matches.get(i);
			if (match.getSpectrum().getId() != previousMatch.getSpectrum().getId()) {
				rank = 1;
			}
			if (match.getScore() == previousMatch.getScore()) {
				rank = previousMatch.getRank();
			}
			match.setRank(rank);
			rank++;
			previousMatch = match;
		}
		//Setting Score ratios for those with rank 1
		int i = matches.size() - 1;
		double previousScore = match.getScore();
		for (; i >= 0; i--) {
			match = matches.get(i);
			if (match.getRank() == 1) {
				match.setScoreRatio(match.getScore() / previousScore);
			} else {
				previousScore = match.getScore();
			}
		}
	}
	
	/**
	 * finds the number of times a certain amino acid is found for each spectrum 
	 * @param matches
	 */
	public static void assignRepeatedPeptideCount(ArrayList<Match> matches) {
		//first error check
		if (matches.size() == 0) return;
		Match.setSortParameter(Match.SORT_BY_SPECTRUM_ID_THEN_PEPTIDE);
		Collections.sort(matches);
		Match match;
		Match previousMatch = matches.get(0);
		int rankCount = 1;
		for (int i = 1; i < matches.size(); i++) {
			//see if these are matches for a different spectrum
			match = matches.get(i);
			if (match.getSpectrum().getId() != previousMatch.getSpectrum().getId()) {
				for (int j = i - rankCount; j < i; j++) {
					matches.get(j).setRepeatCount(rankCount);
				}
				rankCount = 1;
			} else {
				if (match.getPeptide().equals(previousMatch.getPeptide())) {
					rankCount++;
				} else {
					for (int j = i - rankCount; j < i; j++) {
						matches.get(j).setRepeatCount(rankCount);
					}
					rankCount = 1;
				}
			}
			previousMatch = match;
		}
	}
	
	public static void removeDuplicateMatches(ArrayList<Match> matches) {
		//first error check
		if (matches.size() == 0) return;
		
		Match.setSortParameter(Match.SORT_BY_SPECTRUM_ID_THEN_PEPTIDE);
		Collections.sort(matches);
		int numberOfMatches = matches.size();
		Match match;
		Match previousMatch = matches.get(0);
		int spectrumID;
		int previousSpectrumID = previousMatch.getSpectrum().getId();
		boolean areEqual;
		for (int i = 1; i < numberOfMatches; i++) {
			match = matches.get(i);
			spectrumID = match.getSpectrum().getId();
			areEqual = false;
			if (match.getPeptide().isForward() && previousMatch.getPeptide().isForward()) {
				if (match.equals(previousMatch) && match.getPeptide().getStartIndex() == previousMatch.getPeptide().getStartIndex() && spectrumID == previousSpectrumID) {
					areEqual = true;
					matches.remove(i);
					i--;
					numberOfMatches--;
				}
			}
			if (!match.getPeptide().isForward() && !previousMatch.getPeptide().isForward()) {
				if (match.equals(previousMatch) && match.getPeptide().getStartIndex() == previousMatch.getPeptide().getStartIndex() && spectrumID == previousSpectrumID) {
					areEqual = true;
					matches.remove(i);
					i--;
					numberOfMatches--;
				}
			}
			if (!areEqual) {
				previousMatch = match;
				previousSpectrumID = spectrumID;
			}
		}
	}
	
	public static void assignConfidenceValuesToMatches(ArrayList<Match> matches) {
		for (Match match: matches) {
			match.calculateEValue();
			match.calculatePValue();
		}
	}




}
