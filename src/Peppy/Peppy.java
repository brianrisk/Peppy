package Peppy;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;

import HMMScore.HMMScorer;
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
		init(args);
//		splice();
//		runJobs();
		new Peppy(args);
//		cnvThenGenome();
//		cnv();
//		bigJob(args);
//		cnvPeptideMassList();
//		bigHMM(args);
//		testHMMScore();
//		ProteinDigestion.convertDATtoFASTA(Properties.sequenceDirectoryOrFile);
//		runOnProteinDatabase(args);
//		exportPeptideList();
		U.p("done");
	}
	
	public Peppy(String [] args) {
		U.startStopwatch();
		printGreeting();
		peptideTally = 0;
		
		//Load our spectra
		U.p("loading spectra...");
		ArrayList<Spectrum> spectra = Spectrum.loadSpectra();
		U.p("loaded " +spectra.size() + " spectra.");
		
		//Get references to our sequence files -- no nucleotide data is loaded at this point
		ArrayList<Sequence> sequences = Sequence.loadSequences(Properties.sequenceDirectoryOrFile);
		
		//initialize our ArrayList of matches
		ArrayList<Match> matches = new ArrayList<Match>();
		
		if (Properties.useSpliceVariants) {
			//gets the first nucleotide sequence in the first sequence file
			Sequence sequenceFile = sequences.get(0);
			RNA_Sequence rna = new RNA_Sequence(sequenceFile.getNucleotideSequences().get(0), Properties.sequenceRegionStart, Properties.sequenceRegionStop);
			
			U.p("digesting...");
			RNA_Digestor rnaDigestor = new RNA_Digestor(rna);
			ArrayList<Peptide> peptides  = rnaDigestor.getPeptides();
			U.p("peptide tally: " + peptides.size());
			
			U.p("getting matches...");
			matches = getMatches(peptides, spectra, sequenceFile);
			
		} else {
			//loop through each sequence in the sequences ArrayList
			for (Sequence sequence: sequences) {		
				matches.addAll(getMatches(sequence, spectra));
			}
			U.p("peptide tally: " + peptideTally);
		}
	
		U.p("calculating final e values");
		for (Match match: matches) {
			match.calculateEValue();
		}
		
		//create new report directory
		File reportDir = new File(Properties.reportDirectory, Properties.reportDirectoryTitle + " " + System.currentTimeMillis());
		if (Properties.isSequenceFileDNA && Properties.createHTMLReport) {
			U.p("creating HTML reports");
			HTMLReporter report = new HTMLReporter(matches, spectra, sequences, reportDir);
			report.generateFullReport();
		}
		
		U.p("creating text reports");
		TextReporter textReport = new TextReporter(matches, spectra, sequences, reportDir);
		textReport.generateFullReport();
		
		U.p();
		U.stopStopwatch();
	}
	
	public static void runJobs() {
		init();
		File jobsDir = new File("jobs");
		File[] potentialJobsFiles = jobsDir.listFiles();
		ArrayList<File> jobFiles = new ArrayList<File>();
		for (int i = 0; i < potentialJobsFiles.length; i++) {
			if (potentialJobsFiles[i].getName().toLowerCase().endsWith(".txt")) {
				jobFiles.add(potentialJobsFiles[i]);
			}	
		}
		U.p("running " + jobFiles.size() + " jobs");
		for (int i = 0; i < jobFiles.size(); i++) {
			U.p("running job " + (i + 1) + "; " + jobFiles.get(i).getName());
			init();
			Properties.loadProperties(jobFiles.get(i));
			new Peppy(null);
		}
	}

	public static void splice() {
		ArrayList<Sequence> sequences = Sequence.loadSequences(Properties.sequenceDirectoryOrFile);
		for (Sequence sequence: sequences) {
//			RNA_Sequence rna = new RNA_Sequence(sequence.getNucleotideSequences().get(0), 0, sequence.getSequenceLength());
			RNA_Sequence rna = new RNA_Sequence(sequence.getNucleotideSequences().get(0), 35160417, 35253949);
//			rna.printStats();
//			rna.checkCD44Sites();
			U.p("Digesting...");
			RNA_Digestor digestor = new RNA_Digestor(rna);
			U.p(digestor.getPeptides().size());
		}
	}
	
	public static void bigJob(String [] args) {
		String database, error, job;
		for (int i = 0; i < 1; i++) {
			if (i == 0) {
				database = "Genome";
				Properties.isSequenceFileDNA = true;
				Properties.sequenceDirectoryOrFile = new File("/Users/risk2/PeppyOverflow/sequences HG19");
			} else {
				database = "UniProt";
				Properties.isSequenceFileDNA = false;
				Properties.sequenceDirectoryOrFile = new File("/Users/risk2/PeppyOverflow/sequences uniprot/uniprot_trembl_human.dat");
			}
			for (int j = 0; j < 1; j++) {
				if (j == 0) {
					job = "GO_mem_FASP";
					Properties.spectraDirectoryOrFile = new File("/Users/risk2/PeppyOverflow/spectra encode membrane/GO_mem_FASP_dta20100628");
				} else {
					job = "SDS";
					Properties.spectraDirectoryOrFile = new File("/Users/risk2/PeppyOverflow/spectra encode membrane/SDS");
				}
				for (int k = 0; k < 2; k++) {
					if (k == 0) {
						error = "01";
						Properties.spectrumToPeptideMassError = 0.01;
					} else {
						error = "20";
						Properties.spectrumToPeptideMassError = 0.20;
					}
					
					//file:///Volumes/encode/Chris
					Properties.reportDirectory = new File("/Volumes/encode/Chris/report" + "-" + job + "-" + database + "-" + error);
					U.p(Properties.reportDirectory.getName());
					new Peppy(args);
					U.p();
				}
			}
		}
	}
	
	public static void cnvPeptideMassList() {
		//Get references to our sequence files -- no nucleotide data is loaded at this point
		ArrayList<Sequence> sequences = Sequence.loadSequences(Properties.sequenceDirectoryOrFile);
		
		//get sequence, rna, peptides and matches
		Sequence chr11 = sequences.get(0);
		RNA_Sequence rna = new RNA_Sequence(chr11.getNucleotideSequences().get(0), 35160417, 35253949);
		
		U.p("digesting...");
		RNA_Digestor rnaDigestor = new RNA_Digestor(rna);
		ArrayList<Peptide> peptides  = rnaDigestor.getPeptides();
		U.p("peptide tally: " + peptides.size());
		
		File file = new File("CD44-peptides.txt");
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter (new FileWriter(file)));
			int lineTally = 0;
			double previousMass = 0;
			for (Peptide peptide: peptides) {
				if (peptide.getMass() != previousMass) {
					pw.println(peptide.getMass());
					previousMass = peptide.getMass();
					lineTally++;
				}
			}
			U.p("number of unique masses: " + lineTally);
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	/**
	 * This is just to test HMM score
	 */
	public static void testHMMScore() {
		//set up HMM Score
		HMMScore.HMMClass.HmmSetUp();
		
		//Set up peptide, spectra matches;
		ArrayList<Match> matches = new ArrayList<Match>();
		
		matches.add(new Match("/Users/risk2/PeppyOverflow/spectra encode membrane/GO_mem_FASP_dta20100628/020810_M6_FASP_ziptip.dta/020810_M6_FASP_ziptip.9662.9662.2.dta", "ELAEDGYSGVEVR"));
//		matches.add(new Match("HMM-ecoli/CheZ_MSMS_1263.5822_6.txt", "MMDVIQEIER"));
//		matches.add(new Match("HMM-ecoli/CheZ_MSMS_1952.9525_8.txt", "ELGLDQAIAEAAEAIPDAR"));
//		matches.add(new Match("HMM-ecoli/CheZ_MSMS_1642.7903_7.txt", "LYYVVQMTAQAAER"));
//		matches.add(new Match("HMM-ecoli/CheZ_MSMS_1911.8491_5.txt", "ALNSVEASQPHQDQMEK"));
		
		//calculate HMM scores for the spectrum/peptide pairs
		HMMScorer hmmScorer = new HMMScorer(matches);
		hmmScorer.score();
		
		//display the scores
		for (Match match: matches) {
			U.p(match.getPeptide().getAcidSequence() + ": " + match.getScoreHMM());
		}
	}
	
	public static void cnv() {
		U.startStopwatch();
		printGreeting();
		
		U.p("CNV report");
		
		//setting other properties
		Properties.maximumNumberOfMatchesForASpectrum = 2;

		//Load our spectra
		U.p("loading spectra...");
		ArrayList<Spectrum> spectra = Spectrum.loadSpectra();
		U.p("loaded " +spectra.size() + " spectra.");
		
		//Get chr 11
		File chr11File = new File("/Users/risk2/PeppyOverflow/sequences HG19/chr11.fa");
		Sequence chr11 = Sequence.loadSequences(chr11File).get(0);
		RNA_Sequence rna = new RNA_Sequence(chr11.getNucleotideSequences().get(0), 35160417, 35253949);
		
		U.p("digesting...");
		RNA_Digestor rnaDigestor = new RNA_Digestor(rna);
		ArrayList<Peptide> peptides  = rnaDigestor.getPeptides();
		U.p("peptide tally: " + peptides.size());
		
		U.p("getting matches...");
		ArrayList<Match> matches = new ArrayList<Match>();
		matches = getMatches(peptides, spectra, chr11);

		U.p("calculating final e values");
		for (Match match: matches) {
			match.calculateEValue();
		}
		
		File reportDir = new File(Properties.reportDirectory, "report " + System.currentTimeMillis());
		U.p("creating text report");
		TextReporter textReport = new TextReporter(matches, spectra, null, reportDir);
		textReport.generateFullReport();
		
		U.p("creating HTML reports");
		if (Properties.isSequenceFileDNA) {
			HTMLReporter report = new HTMLReporter(matches, spectra, null, reportDir);
			report.generateFullReport();
		}

		U.p();
		U.stopStopwatch();
	}
	
	public static void cnvThenGenome() {
		U.startStopwatch();
		printGreeting();
		
		U.p("Searching the bit of chr11, then comparing to genome");
		
		//setting other properties
		Properties.maximumNumberOfMatchesForASpectrum = 2;

		//Load our spectra
		U.p("loading spectra...");
		ArrayList<Spectrum> spectra = Spectrum.loadSpectra();
		U.p("loaded " +spectra.size() + " spectra.");
		
		//get sequence, rna, peptides and matches
		File chr11File = new File("/Users/risk2/PeppyOverflow/sequences HG19/chr11.fa");
		Sequence chr11 = Sequence.loadSequences(chr11File).get(0);
		RNA_Sequence rna = new RNA_Sequence(chr11.getNucleotideSequences().get(0), 35160417, 35253949);
		
		U.p("digesting...");
		RNA_Digestor rnaDigestor = new RNA_Digestor(rna);
		ArrayList<Peptide> peptides  = rnaDigestor.getPeptides();
		U.p("peptide tally: " + peptides.size());
		
		U.p("getting matches...");
		ArrayList<Match> chr11matches = new ArrayList<Match>();
		chr11matches = getMatches(peptides, spectra, chr11);

		U.p("calculating final e values");
		for (Match match: chr11matches) {
			match.calculateEValue();
		}
		
		//clearing old peptides
		peptides = null;
		System.gc();
		
		//finding reduced set of promising spectra
		ArrayList<Spectrum> interestingSpectra = new ArrayList<Spectrum>();
		for (Match match: chr11matches) {
			if (match.getRank() == 1 && match.getEValue() > Properties.eValueCutOff) {
				match.getSpectrum().clearEValues();
				interestingSpectra.add(match.getSpectrum());
			}
		}
		
		//Search genome
		ArrayList<Sequence> sequences = Sequence.loadSequences(new File("/Users/risk2/PeppyOverflow/sequences HG19"));
		ArrayList<Match> genomeMatches = new ArrayList<Match>();
		for (Sequence sequence: sequences) {		
			genomeMatches.addAll(getMatches(sequence, interestingSpectra));
		}
		
		//find matches from chr11 that outmatch the whole genome
		ArrayList<Match> interestingMatches = new ArrayList<Match>();
		for (Match match: chr11matches) {
			for (Match genomeMatch: genomeMatches) {
				if (match.getSpectrum().getId() == genomeMatch.getSpectrum().getId()) {
					if (match.getScore() >= genomeMatch.getScore()) {
						interestingMatches.add(match);
					}
				}
			}
		}
		File reportDir = new File(Properties.reportDirectory, "report " + System.currentTimeMillis());
		
		U.p("creating text report");
		TextReporter textReport = new TextReporter(interestingMatches, interestingSpectra, sequences, reportDir);
		textReport.generateFullReport();
		
		U.p("creating HTML reports");
		if (Properties.isSequenceFileDNA) {
			HTMLReporter report = new HTMLReporter(interestingMatches, interestingSpectra, sequences, reportDir);
			report.generateFullReport();
		}

		U.p();
		U.stopStopwatch();
	}
	

	public static void init(String propertiesFile) {
		System.setProperty("java.awt.headless", "true"); 
		Properties.loadProperties(propertiesFile);
		if (Properties.defaultScore == Properties.DEFAULT_SCORE_HMM) {
			HMMScore.HMMClass.HmmSetUp();
		}
	}
	
	public static void init(String [] args) {
		if (args.length == 0) {
			init();
		} else {
			init(args[0]);
		}
	}
	
	public static void init() {
		init("properties.txt");
	}
	
	
	public static ArrayList<Match> getMatches(Sequence sequence, ArrayList<Spectrum> spectra) {
			U.p("working on sequence: " +sequence.getSequenceFile().getName());
			ArrayList<Match> matches = new ArrayList<Match>() ;
			ArrayList<Peptide> peptides;
			if (Properties.isSequenceFileDNA) {
				peptides = sequence.extractMorePeptides();
				//continually extract peptides from the sequence until there aren't anymore
				while (peptides != null) {
					peptideTally += peptides.size();
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
				peptideTally += peptides.size();
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
			assignRankCountsToMatches(matches);

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
		peptideTally += peptides.size();
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
	}
	
	public static void assignRankCountsToMatches(ArrayList<Match> matches) {
		//first error check
		if (matches.size() == 0) return;
		Match.setSortParameter(Match.SORT_BY_SPECTRUM_ID_THEN_SCORE);
		Collections.sort(matches);
		Match match;
		Match previousMatch = matches.get(0);
		int rankCount = 1;
		for (int i = 1; i < matches.size(); i++) {
			//see if these are matches for a different spectrum
			match = matches.get(i);
			if (match.getSpectrum().getId() != previousMatch.getSpectrum().getId()) {
				for (int j = i - rankCount; j < i; j++) {
					matches.get(j).setRankCount(rankCount);
				}
				rankCount = 1;
			} else {
				if (match.getScore() == previousMatch.getScore()) {
					rankCount++;
				} else {
					for (int j = i - rankCount; j < i; j++) {
						matches.get(j).setRankCount(rankCount);
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
