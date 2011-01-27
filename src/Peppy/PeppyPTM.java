package Peppy;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;

import Reports.TextReporter;
import Utilities.U;




/**
 * Peppy
 * A very stripped down Java version of the MS/MS to genome mapping of Morgan Gidding's GFS.
 * Designed with the following goals:
 * 1) More simple code to promote open source development
 * 2) Takes advantage of Java's popularity over Objective-C
 * 3) better multi-threading
 * @author Brian Risk
 *
 */
public class PeppyPTM extends Peppy{
	
	
	public static void main(String [] args) {
		init(args);
		new PeppyPTM(args);
		U.p("done");
	}
	
	public PeppyPTM(String [] args) {
		super(args);
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
		for (Sequence sequence: sequences) {
			proteins.addAll(ProteinDigestion.getProteinsFromDatabase(sequence.getSequenceFile()));
		}
		
		//digest peptides from all those proteins
		peptides = ProteinDigestion.getPeptidesFromListOfProteins(proteins);
		
		//get the matches
		matches = getMatches(peptides, spectra, null);
		
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
		
		//get the peptides that weren't found
		ArrayList<Peptide> unfoundPeptides = new ArrayList<Peptide>();
		for (Protein protein: topProteins) {
			unfoundPeptides.addAll(protein.getUnfoundPeptides());
		}
		
		//create new report directory
		File reportDir = new File(Properties.reportDirectory, Properties.reportDirectoryTitle + " " + System.currentTimeMillis());
		
		U.p("creating text reports");
		TextReporter textReport = new TextReporter(matches, spectra, sequences, reportDir);
		textReport.generateFullReport();
		textReport.generatePropertiesFile();
		
		U.p();
		U.stopStopwatch();
	}
	
	
	

}
