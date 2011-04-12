package Peppy;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;

import Reports.PeptidesWithModificaitonsHTMLPage;
import Reports.ProteinsHTMLPage;
import Utilities.U;


/**
 * PeppyPTM
 * An extension of Peppy to include protein modifications
 * @author Brian Risk
 *
 */
public class Peppy_VariMod extends Peppy{
	
	/*
	public static void main(String [] args) {
		init(args);
		runPeppyPTM(args);
		U.p("done");
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
		ArrayList<Sequence> sequences = Sequence.loadSequenceFiles(Properties.sequenceDirectoryOrFile);
		
		//initialize our ArrayList of matches
		ArrayList<Match> unModifiedMatches = new ArrayList<Match>();
		
		//Set up our proteins
		ArrayList<Protein> proteins = new ArrayList<Protein>();
		
		//Set up our peptides
		ArrayList<Peptide> proteinPeptides;
		
		//loop through sequences, get proteins
		U.p("loading proteins");
		for (Sequence sequence: sequences) {
			proteins.addAll(Sequence_Protein.getProteinsFromDatabase(sequence.getSequenceFile()));
		}
		
		//digest peptides from all those proteins
		U.p("digesting proteins");
		proteinPeptides = Sequence_Protein.getPeptidesFromListOfProteins(proteins);
		U.p(proteinPeptides.size() + " peptides created");
		
		//get the matches
		U.p("finding matches");
		unModifiedMatches = getMatchesWithPeptides(proteinPeptides, spectra);
		
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
		int topProteinNumber = 64;
		if (topProteinNumber > proteins.size()) topProteinNumber = proteins.size();
		for (int i = 0; i < topProteinNumber; i++) {
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
		ArrayList<Match_IMP_VariMod> matchesPTM = (new ScoringThreadServerPTM(unfoundPeptides, spectra)).getMatches();
		Collections.sort(matchesPTM);
		
		//find matches with greater than 2 Da difference, and great score
		//TODO why are some from less than 2 Da not showing up in first search?
		Match_IMP_VariMod topMatch;
		ArrayList<Match_IMP_VariMod> topMatchesPTM = new ArrayList<Match_IMP_VariMod>();
		for (int i = 0; i < matchesPTM.size(); i++) {
			topMatch = matchesPTM.get(i);
			if (topMatch.getModificationMass() > 2) {
				if (topMatch.getScore() >= 21) {
					topMatchesPTM.add(topMatch);
				}
			}
		}
		
		
		//add the good matches we've found to our proteins
		for (Match_IMP_VariMod match: topMatchesPTM) {
			match.getPeptide().getProtein().addMatchPTM(match);
		}
		
		U.p("creating reports");
		//create report directory
		File reportDir = new File(Properties.reportDirectory, "Proteins " + System.currentTimeMillis());
		reportDir.mkdirs();
		
		//making the main index page
		ProteinsHTMLPage php = new ProteinsHTMLPage(topProteins, new File(reportDir, "index.html"));
		php.makePage();
		
		//the modifications page
		PeptidesWithModificaitonsHTMLPage pwmhp = new PeptidesWithModificaitonsHTMLPage(topProteins, new File(reportDir, "modifications.html"));
		pwmhp.makePage();
		
		
		U.p();
		U.stopStopwatch();
	}
	*/
	
	

}
