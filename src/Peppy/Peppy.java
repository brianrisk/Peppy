package Peppy;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;

import Reports.HTMLReporter;
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
public class Peppy {
	
	//So that we may report the total amount of peptides found
	static int peptideTally = 0;
	
	public static void main(String [] args) {
		init();
//		new Peppy(args);
		testHMMScore();
//		ProteinDigestion.convertDATtoFASTA(Properties.sequenceDirectoryOrFile);
//		runOnProteinDatabase(args);
//		exportPeptideList();
		U.p("done");
	}
	

	
	/**
	 * This is just to test HMM score
	 */
	public static void testHMMScore() {
		//set up HMM Score
		HMMScore.HMMClass.HmmSetUp();
		
		//Set up peptide, spectra matches;
		ArrayList<Match> matches = new ArrayList<Match>();
		matches.add(new Match("HMM-ecoli/CheZ_MSMS_1263.5822_6.txt", "MMDVIQEIER"));
		matches.add(new Match("HMM-ecoli/CheZ_MSMS_1952.9525_8.txt", "ELGLDQAIAEAAEAIPDAR"));
		matches.add(new Match("HMM-ecoli/CheZ_MSMS_1642.7903_7.txt", "LYYVVQMTAQAAER"));
		matches.add(new Match("HMM-ecoli/CheZ_MSMS_1911.8491_5.txt", "ALNSVEASQPHQDQMEK"));
		
		//calculate HMM scores for the spectrum/peptide pairs
		HMMScorer hmmScorer = new HMMScorer(matches);
		hmmScorer.score();
		
		//display the scores
		for (Match match: matches) {
			U.p(match.getPeptide().getAcidSequence() + ": " + match.getScoreHMM());
		}
	}

	public Peppy(String [] args) {
			U.startStopwatch();
			printGreeting();
			
			
			//setting other properties
			Properties.maximumNumberOfMatchesForASpectrum = 2;
			Properties.reduceDuplicateMatches = false;
	//		Properties.peakDifferenceThreshold = 0.25;
	//		Properties.peakIntensityExponent = 0.3;
//			Properties.spectrumToPeptideMassError = 0.01;
			Properties.spectrumToPeptideMassError = 0.1;
			Properties.peakDifferenceThreshold = 0.8;
			
			//Load our spectra
			U.p("loading spectra...");
			ArrayList<Spectrum> spectra = Spectrum.loadSpectra();
			U.p("loaded " +spectra.size() + " spectra.");
			
			//Get references to our sequence files -- no nucleotide data is loaded at this point
			ArrayList<Sequence> sequences = Sequence.loadSequences(Properties.sequenceDirectoryOrFile);
			
			//initialize our ArrayList of matches
			ArrayList<Match> matches = new ArrayList<Match>();
			
			//loop through each sequence in the sequences ArrayList
			for (Sequence sequence: sequences) {		
				matches.addAll(getMatches(sequence, spectra));
			}
			U.p("peptide tally: " + peptideTally);

			//recalculate e values now that all peptides are in for the count
			U.p("calculating final e values");
			for (Match match: matches) {
				match.calculateEValue();
			}
			
			ArrayList<Match> specificMatches =  new ArrayList<Match>();
			for (Match match: matches) {
				int index = match.getPeptide().getIndex();
				if (index >= 31862431 && index <= 38471004) {
					specificMatches.add(match);
				}
			}
			matches = specificMatches;
			
			//calculate HMM scores
//			HMMScorer hmmScorer = new HMMScorer(matches);
//			hmmScorer.score();
			
			U.p("Creating reports");
//			if (Properties.isSequenceFileDNA) {
				HTMLReporter report = new HTMLReporter(matches, spectra, sequences);
				report.generateFullReport();
//			}
			
			TextReporter textReport = new TextReporter(matches, spectra, sequences);
			textReport.generateFullReport();
			
			U.p();
			U.stopStopwatch();
		}

	public static void init() {
		System.setProperty("java.awt.headless", "true"); 
		Properties.loadProperties("properties.txt");
		HMMScore.HMMClass.HmmSetUp();
	}
	
	/**
	 * 
	 * @param args an array of arguments.  Element 0 is the spectral directory; Element 1 is the FASTA database.
	 */
	public static void runOnProteinDatabase(String [] args) {
		U.startStopwatch();
		printGreeting();	
		
		if (args[0] == null || args[0].equals("")) {
			U.p();
			U.p("Here is how to use it:");
			U.p("All file paths can be absolute or relative to where this app is located.");
			U.p("Arg 0: The file path directory to the folder which contains the spectra.");
			U.p("Arg 1: The file path directory to the FASTA protein database file.");
		} else {
			Properties.spectraDirectoryOrFile = new File(args[0]);
			
			//setting properties
			Properties.maximumNumberOfMatchesForASpectrum = 2;
			Properties.peakDifferenceThreshold = 0.25;
			Properties.peakIntensityExponent = 0.3;
			
			//Load our spectra
			ArrayList<Spectrum> spectra = Spectrum.loadSpectra();
			U.p("loaded " +spectra.size() + " spectra.");
	
	
			//load the peptides from the database
			ArrayList<Peptide> peptides = ProteinDigestion.getPeptidesFromFASTA(new File(args[1]));
			
			
			ArrayList<Match> matches = (new ScoringEngine(peptides, spectra, null)).getMatches();
			
			TextReporter textReport = new TextReporter(matches, spectra, null);
			textReport.generateFullReport();
			
			U.p();
			U.stopStopwatch();
		}
	}
	
	
	public static ArrayList<Match> getMatches(Sequence sequence, ArrayList<Spectrum> spectra) {
			U.p("Working on sequence: " +sequence.getSequenceFile().getName());
			ArrayList<Match> matches = new ArrayList<Match>() ;
			ArrayList<Peptide> peptides;
			if (Properties.isSequenceFileDNA) {
				peptides = sequence.extractMorePeptides();
				//continually extract peptides from the sequence until there aren't anymore
				while (peptides != null) {
					peptideTally += peptides.size();
					//This is where the bulk of the processing in long jobs takes
					ArrayList<Match> newMatches = (new ScoringEngine(peptides, spectra, sequence)).getMatches();
					//Add only matches with a decent e value
					for (Match match: newMatches) {
						if (match.getEValue() <= Properties.eValueCutOff) matches.add(match);
					}
//					matches.addAll(newMatches);
					//free up the memory of the old peptide arraylist
					peptides.clear();
					System.gc();
					peptides = sequence.extractMorePeptides();
				}
				sequence.clearNucleotideData();
			} else {
				peptides = ProteinDigestion.getPeptidesFromDatabase(sequence.getSequenceFile());
				peptideTally += peptides.size();
				//This is where the bulk of the processing in long jobs takes
				ArrayList<Match> newMatches = (new ScoringEngine(peptides, spectra, sequence)).getMatches();
				//Add only matches with a decent e value
				for (Match match: newMatches) {
					if (match.getEValue() <= Properties.eValueCutOff) matches.add(match);
				}
//				matches.addAll(newMatches);
			}
			return matches;
	}
	
	
	/**
	 * saves a file of peptide information
	 * @param peptides
	 */
	@SuppressWarnings("unused")
	private static void exportPeptideList() {
		//Get references to our sequence files -- no nucleotide data is loaded at this point
		ArrayList<Sequence> sequences = Sequence.loadSequences(Properties.sequenceDirectoryOrFile);
		
		try {
			//loop through each sequence in the sequences ArrayList
			File peptidesFolder = new File ("peptides");
			peptidesFolder.mkdir();
			for (int sequenceIndex = 0; sequenceIndex < sequences.size(); sequenceIndex++) {
				Sequence sequence = sequences.get(sequenceIndex);
				PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(peptidesFolder, sequence.getSequenceFile().getName()))));
				
				ArrayList<Peptide> peptides = null;
				if (Properties.isSequenceFileDNA) {
					//
				} else {
					peptides = ProteinDigestion.getPeptidesFromDatabase(sequence.getSequenceFile());
				}
				
				for (int i = 1; i < peptides.size(); i++) {
					pw.println(peptides.get(i));
				}
				peptides = null;
				System.gc();
				pw.flush();
				pw.close();
				U.p(sequence.getSequenceFile().getName() + " digested.");
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private static void printGreeting() {
		U.p("Welcome to Peppy");
		U.p("Proteogenomic mapping software.");
		U.p("Developed 2010 by the Giddings Lab");
		U.p();
	}
	

}
