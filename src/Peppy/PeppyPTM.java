package Peppy;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

import Reports.ProteinReporter;
import Reports.ProteinsHTMLPage;
import SpectralVisualizer.SpectralVisualizerPTM;
import Utilities.U;




/**
 * PeppyPTM
 * An extension of Peppy to include protein modifications
 * @author Brian Risk
 *
 */
public class PeppyPTM extends Peppy{
	
	
	public static void main(String [] args) {
		init(args);
//		test();
		runPeppyPTM(args);
//		runPeppyProtein(args);
//		findPeptidePTM();
		U.p("done");
	}
	
//	public static void test() {
//		ArrayList<Spectrum> spectra = Spectrum.loadSpectra();
//		Spectrum spectrum = spectra.get(0);
//		U.p(spectrum.getFile().getName());
//		Peptide peptide = new Peptide("PSYVLSGAAMNVAYTDGDLER");
//		MatchPTM matchPTM = new MatchPTM(spectrum, peptide);
//		U.p("score: " + matchPTM.getScore());
//		U.p("difference: " + matchPTM.difference);
//		
//		File visualizationTestFile = new File("visualizationTest.jpg");
//		try {
//			SpectralVisualizer.drawDeluxSpectrum(spectrum, peptide, visualizationTestFile);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//	}
	
	public static void findPeptidePTM() {
		U.p("Peppy+ finding PTM");
		U.startStopwatch();
		
		//Load our spectra
		U.p("loading spectra...");
		ArrayList<Spectrum> spectra = Spectrum.loadSpectra();
		U.p("loaded " +spectra.size() + " spectra.");
		Collections.sort(spectra);
		

		//initialize our ArrayList of matches
		ArrayList<MatchPTM> matches = new ArrayList<MatchPTM>();
		
		
		//Set up our peptides
		ArrayList<Peptide> peptides = new ArrayList<Peptide>();
//		Peptide thePeptide = new Peptide("AAHSEGNTTAGLDMR");
		Peptide thePeptide = new Peptide("LYVGNIPFGITEEAMMDFFNAQMR");
		peptides.add(thePeptide);
		
		matches = (new ScoringThreadServerPTM(peptides, spectra)).getMatches();
		Collections.sort(matches);
		U.p("number of matches: " + matches.size());
		
		MatchPTM match = matches.get(0);
		Spectrum spectrum = match.getSpectrum();
		Peptide peptide = match.getPeptide();
		String acidString = peptide.getAcidSequenceString();
		U.p("higest matching spectum: " +match.getSpectrum().getFile().getName());
		
		double difference = ( match.getSpectrum().getMass() - match.getPeptide().getMass());
		U.p("difference: " + difference);
		U.p("raw score: " + match.getScore());
		U.p("second score: " + matches.get(1).getScore());
		
		double imp;
		double bestIMP = Double.MAX_VALUE;
		int bestIndex = 0;
		MatchPTM matchPTM;
		for (int i= 0; i < acidString.length(); i++) {
			matchPTM = new MatchPTM(spectrum, peptide);
			imp = matchPTM.calculateIMP(matchPTM.difference, i);
			U.p(i + " " + acidString.charAt(i) + ": " + imp);
			if (imp < bestIMP) {
				bestIMP = imp;
				bestIndex = i;
			}
		}
		
		try {
			SpectralVisualizerPTM.drawDeluxSpectrum(spectrum, peptide, new File ("PTM_IS_AWESOME.jpg"), difference, bestIndex);
		} catch (IOException e) {
			e.printStackTrace();
		}
		U.p();
		U.stopStopwatch();
	}
	
	
	
	public static void runPeppyProtein(String [] args) {
		U.p("Peppy+ standard protein report");
		U.startStopwatch();
		
		//Load our spectra
		U.p("loading spectra...");
		ArrayList<Spectrum> spectra = Spectrum.loadSpectra();
		U.p("loaded " +spectra.size() + " spectra.");
		
		//Get references to our sequence files
		//These should be protein databases such as UniProt
		ArrayList<Sequence> sequences = Sequence.loadSequences(Properties.sequenceDirectoryOrFile);
		
		//initialize our ArrayList of matches
		ArrayList<Match> matches = new ArrayList<Match>();
		
		//Set up our proteins
		ArrayList<Protein> proteins = new ArrayList<Protein>();
		
		//Set up our peptides
		ArrayList<Peptide> peptides;
		
		//loop through sequences, get proteins
		U.p("loading proteins");
		for (Sequence sequence: sequences) {
			proteins.addAll(ProteinDigestion.getProteinsFromDatabase(sequence.getSequenceFile()));
		}
		
		//digest peptides from all those proteins
		U.p("digesting proteins");
		peptides = ProteinDigestion.getPeptidesFromListOfProteins(proteins);
		
		//get the matches
		U.p("finding matches");
		matches = getMatches(peptides, spectra, null);
		
		//double check e values
		assignConfidenceValuesToMatches(matches);
		
		//add the good matches we've found to our proteins
		for (Match match: matches) {
			match.getPeptide().getProtein().addMatch(match);
		}
		
		//Sorting our proteins by the score they have now acquired
		Collections.sort(proteins);
		
		//get the top proteins
		ArrayList<Protein> topProteins = new ArrayList<Protein>();
		for (int i = 0; i < 64; i++) {
			topProteins.add(proteins.get(i));
		}
		
//		//get the peptides that weren't found
//		ArrayList<Peptide> unfoundPeptides = new ArrayList<Peptide>();
//		for (Protein protein: topProteins) {
//			unfoundPeptides.addAll(protein.getUnfoundPeptides());
//		}
//		
//		//quick update on things
//		U.p("there were this many unfound peptides: " + unfoundPeptides.size());
		
		//create new report directory
		File reportDir = new File(Properties.reportDirectory, "Proteins " + System.currentTimeMillis());
		reportDir.mkdirs();
		ProteinReporter reporter = new ProteinReporter(topProteins, new File(reportDir, "index.html"));
		reporter.generateReport();
		
		U.beep();
		U.p();
		U.stopStopwatch();
	}
	
	public static void runPeppyPTM(String [] args) {
		U.p("Peppy+ PTM");
		U.startStopwatch();
		
		//Load our spectra
		U.p("loading spectra...");
		ArrayList<Spectrum> spectra = Spectrum.loadSpectra();
		U.p("loaded " +spectra.size() + " spectra.");
		
		//Get references to our sequence files
		//These should be protein databases such as UniProt
		ArrayList<Sequence> sequences = Sequence.loadSequences(Properties.sequenceDirectoryOrFile);
		
		//initialize our ArrayList of matches
		ArrayList<Match> unModifiedMatches = new ArrayList<Match>();
		
		//Set up our proteins
		ArrayList<Protein> proteins = new ArrayList<Protein>();
		
		//Set up our peptides
		ArrayList<Peptide> proteinPeptides;
		
		//loop through sequences, get proteins
		U.p("loading proteins");
		for (Sequence sequence: sequences) {
			proteins.addAll(ProteinDigestion.getProteinsFromDatabase(sequence.getSequenceFile()));
		}
		
		//digest peptides from all those proteins
		U.p("digesting proteins");
		proteinPeptides = ProteinDigestion.getPeptidesFromListOfProteins(proteins);
		
		//get the matches
		U.p("finding matches");
		unModifiedMatches = getMatches(proteinPeptides, spectra, null);
		
		//double check e values
		assignConfidenceValuesToMatches(unModifiedMatches);
		
		//add the good matches we've found to our proteins
		for (Match match: unModifiedMatches) {
			match.getPeptide().getProtein().addMatch(match);
		}
		
		//Sorting our proteins by the score they have now acquired
		Collections.sort(proteins);
		
		//get the top proteins
		ArrayList<Protein> topProteins = new ArrayList<Protein>();
		for (int i = 0; i < 64; i++) {
			topProteins.add(proteins.get(i));
		}
		
		//get the peptides that weren't found
		ArrayList<Peptide> unfoundPeptides = new ArrayList<Peptide>();
		for (Protein protein: topProteins) {
			unfoundPeptides.addAll(protein.getUnfoundPeptides());
		}
		
		//quick update on things
		U.p("there were this many unfound peptides: " + unfoundPeptides.size());
		
		//search for all unfound peptides
		U.p("finding modifications");
		Properties.maximumNumberOfMatchesForASpectrum = 1;
		ArrayList<MatchPTM> matchesPTM = (new ScoringThreadServerPTM(unfoundPeptides, spectra)).getMatches();
		Collections.sort(matchesPTM);
		
		//find matches with greater than 2 Da difference, and great score
		//TODO why are some from less than 2 Da not showing up in first search?
		MatchPTM topMatch;
		ArrayList<MatchPTM> topMatchesPTM = new ArrayList<MatchPTM>();
		for (int i = 0; i < matchesPTM.size(); i++) {
			topMatch = matchesPTM.get(i);
			if (topMatch.getDifference() > 2) {
				if (topMatch.getScore() >= 21) {
					topMatchesPTM.add(topMatch);
				}
			}
		}
		
		
		//add the good matches we've found to our proteins
		for (MatchPTM match: topMatchesPTM) {
			match.getPeptide().getProtein().addMatchPTM(match);
		}
		
		//create report
		U.p("creating reports");
		File reportDir = new File(Properties.reportDirectory, "Proteins " + System.currentTimeMillis());
		reportDir.mkdirs();
		ProteinsHTMLPage php = new ProteinsHTMLPage(topProteins, new File(reportDir, "index.html"));
		php.makePage();
		
		
		U.p();
		U.stopStopwatch();
	}
	
	
	

}
