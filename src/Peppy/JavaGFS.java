package Peppy;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;

import Reports.TextReporter;
import Utilities.U;




/**
 * JavaGFS
 * A very stripped down Java version of Morgan Gidding's GFS.
 * Designed with the following goals:
 * 1) More simple code to promote open source development
 * 2) Takes advantage of Java's popularity over Objective-C
 * 3) better multi-threading
 * @author Brian Risk
 *
 */
public class JavaGFS {
	
	
	public static void main(String [] args) {
		System.setProperty("java.awt.headless", "true"); 
		Properties.loadProperties("properties.txt");
//		ReliabilityTester.runTheGamut("USP");
//		ReliabilityTester.runTheGamut("ecoli");
//		ReliabilityTester.testReliability("aurum");
//		ReliabilityTester.testReliability("human");
//		ReliabilityTester.testReliability("ecoli");
//		ReliabilityTester.exportHighScoringPeptidesFromSwissProt("USP");
//		ReliabilityTester.exportPeptidesInCommonWithDatabase("human");
//		ReliabilityTester.exportPeptidesInCommonWithDatabase("ecoli");
//		ReliabilityTester.makeSureWeAreProperlyDigestingTheGenome("ecoli");
		new JavaGFS	(args);
//		testHMMScoreOnSwissProt();
//		runOnProteinDatabase(args);
//		exportPeptideList();
		U.p("done");
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
			Properties.spectraFile = new File(args[0]);
			
			//setting properties
			Properties.maximumNumberOfMatchesForASpectrum = 1;
			Properties.peakDifferenceThreshold = 0.25;
			Properties.peakIntensityExponent = 0.3;
			
			//Load our spectra
			ArrayList<Spectrum> spectra = Spectrum.loadSpectra();
			U.p("loaded " +spectra.size() + " spectra.");
	
	
			//load the peptides from the database
			ArrayList<Peptide> peptides = ProteinDigestion.getPeptidesFromProteinFile(new File(args[1]));
			
			
			ArrayList<SpectrumPeptideMatch> matches = JavaGFS.asynchronousDigestion(peptides, spectra, null);
			
			TextReporter textReport = new TextReporter(matches, spectra, null);
			textReport.generateFullReport();
			
			U.p();
			U.stopStopwatch();
		}
	}
	
	
	/**
	 * This is just to test HMM score on swiss prot
	 */
	public static void testHMMScoreOnSwissProt() {
		U.startStopwatch();
		printGreeting();	

//		Properties.spectraFile = new File("spectra human677");
		Properties.spectraFile = new File("spectra USP");
//		Properties.spectraFile = new File("spectra problem");
		
		//setting properties
		Properties.maximumNumberOfMatchesForASpectrum = 1;
		Properties.peakDifferenceThreshold = 0.25;
		Properties.peakIntensityExponent = 0.3;
		
		//Load our spectra
		ArrayList<Spectrum> spectra = Spectrum.loadSpectra();
		U.p("loaded " +spectra.size() + " spectra.");

		//set up HMM
		HMMScore.HMMClass.HmmSetUp();

		//load the peptides from the database
//		ArrayList<Peptide> peptides = ProteinDigestion.getPeptidesFromProteinFile(new File("tests/databases/uniprot_sprot.fasta"));
		ArrayList<Peptide> peptides = ProteinDigestion.getPeptidesFromProteinFile(new File("tests/databases/ipi.HUMAN.v3.53.fasta"));
		
		
		ArrayList<SpectrumPeptideMatch> matches = JavaGFS.asynchronousDigestion(peptides, spectra, null);
		
		TextReporter textReport = new TextReporter(matches, spectra, null);
		textReport.generateFullReport();
		
		U.p();
		U.stopStopwatch();
	}
	
	
	
	public JavaGFS(String [] args) {
		U.startStopwatch();
		printGreeting();
		
		
		//setting other properties
//		Properties.maximumNumberOfMatchesForASpectrum = 1;
//		Properties.peakDifferenceThreshold = 0.25;
//		Properties.peakIntensityExponent = 0.3;
		
		//Load our spectra
		ArrayList<Spectrum> spectra = Spectrum.loadSpectra();
		U.p("loaded " +spectra.size() + " spectra.");
		
		//Get references to our sequence files -- no nucleotide data is loaded at this point
		ArrayList<Sequence> sequences = Sequence.loadSequences(Properties.sequenceFile);
		
		//initialize our ArrayList of matches
		ArrayList<SpectrumPeptideMatch> matches = new ArrayList<SpectrumPeptideMatch>();
		
		//set up HMM
//		HMMScore.HMMClass.HmmSetUp();
		
		//loop through each sequence in the sequences ArrayList
		//for (int sequenceIndex = 0; sequenceIndex < sequences.size(); sequenceIndex++) {
		for (Sequence sequence: sequences) {		
			matches.addAll(asynchronousDigestion(sequence, spectra));
		}
		
		//calculate HMM scores
//		HMMScorer hmmScorer = new HMMScorer(matches);
//		hmmScorer.score();
		
//		HTMLReporter report = new HTMLReporter(matches, spectra, sequences);
//		report.generateFullReport();
		
		TextReporter textReport = new TextReporter(matches, spectra, sequences);
		textReport.generateFullReport();
		
//		saveBadSpectra(matches);
		
		U.p();
		U.stopStopwatch();
	}
	
	public static void saveBadSpectra(ArrayList<SpectrumPeptideMatch> matches) {
		try {
			File badSpectraFolder = new File("bad spectra");
			badSpectraFolder.mkdirs();
			File badPeptidesFile = new File("bad-peptides.txt");
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(badPeptidesFile)));
			int fileIndex = 1;
			for (int i = 0; i < matches.size(); i++) {
				SpectrumPeptideMatch match = matches.get(i);
				if (match.getScoreHMM() < 1) {
					//copy spectrum file
					InputStream in = new FileInputStream(match.getSpectrum().getFile());
					File badSpectrum = new File(badSpectraFolder, fileIndex + ".dta");	
					OutputStream out = new FileOutputStream(badSpectrum);
					byte[] buf = new byte[1024];
					int len;
					while ((len = in.read(buf)) > 0) {
						out.write(buf, 0, len);
					}
					in.close();
					out.close();
					fileIndex++;
					
					//save bad peptide
					pw.println(match.getPeptide().getAcidSequence());
				}
				
			}
			pw.flush();
			pw.close();
		}
		catch (FileNotFoundException fnfe) {
			fnfe.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/**
	 * For the sake of a lower memory footprint we will not extract
	 * every reading frame, every missed and non-missed cleavages combination
	 * all at once.  For decent sized chromosomes this would require many gigabytes
	 * of memory.  Many users will not have enough.  To that end we will process
	 * the genome in pieces.  One reading frame at a time.  Once with missed cleavages
	 * and once without.
	 */
	public static ArrayList<SpectrumPeptideMatch> synchronousDigestion(Sequence sequence, ArrayList<Spectrum> spectra) {
		ArrayList<SpectrumPeptideMatch> matches = new ArrayList<SpectrumPeptideMatch>();
		ArrayList<NucleotideSequence> nucleotideSequences = sequence.getNucleotideSequences();
		
		//ArrayList<Peptide> peptides = sequence.extractPeptides();
		for (byte frame = 0; frame < 3; frame++) {
			for (int forwards = 0; forwards < 2; forwards++) {
				for (int missedCleavage = 0; missedCleavage < 2; missedCleavage++) {
					//we extract our list of peptides
					ArrayList<Peptide> peptides = sequence.extractPeptides(nucleotideSequences, frame, forwards == 0, missedCleavage == 0);
					
					//This is where the bulk of the processing in long jobs takes
					ScoringEngine engine = new ScoringEngine(peptides, spectra, sequence);
					
					//harvest the results
					matches.addAll(engine.getMatches());
					
					//I know Java doesn't need memory management and all that, but this is a lot of memory we're talkin' about here
					peptides = null;
					System.gc();
				}
			}
		}
		return matches;
	}
	
	/**
	 * takes full advantage of the SequenceDigestionThread.  However, this 
	 * method requires much more memory.
	 */
	public static ArrayList<SpectrumPeptideMatch> asynchronousDigestion(Sequence sequence, ArrayList<Spectrum> spectra) {
			//This is where the big memory drain comes from.  We are extracting
			//a list of peptides from the sequence file.
			U.p("Digesting file: " +sequence.getSequenceFile().getName());
			ArrayList<Peptide> peptides;
			if (Properties.isSequenceFileDNA) {
				peptides = sequence.extractPeptides();
			} else {
				peptides = ProteinDigestion.getPeptidesFromProteinFile(sequence.getSequenceFile());
			}
			
			return asynchronousDigestion(peptides, spectra, sequence);
	}
	
	/**
	 * takes full advantage of the SequenceDigestionThread.  However, this 
	 * method requires much more memory.
	 */
	public static ArrayList<SpectrumPeptideMatch> asynchronousDigestion(ArrayList<Peptide> peptides, ArrayList<Spectrum> spectra, Sequence sequence) {

			//This is where the bulk of the processing in long jobs takes
			ScoringEngine engine = new ScoringEngine(peptides, spectra, sequence);
			
			//harvest the results
			ArrayList<SpectrumPeptideMatch> matches =  engine.getMatches();
			
			//I know Java doesn't need memory management and all that, but let's just be sure
			peptides = null;
			System.gc();
			
			return matches;
	}
	
	
	/**
	 * saves a file of peptide information
	 * @param peptides
	 */
	@SuppressWarnings("unused")
	private static void exportPeptideList() {
		//Get references to our sequence files -- no nucleotide data is loaded at this point
		ArrayList<Sequence> sequences = Sequence.loadSequences(Properties.sequenceFile);
		
		try {
			//loop through each sequence in the sequences ArrayList
			for (int sequenceIndex = 0; sequenceIndex < sequences.size(); sequenceIndex++) {
				Sequence sequence = sequences.get(sequenceIndex);
				PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("peptides/" + sequence.getSequenceFile().getName())));
				
				ArrayList<Peptide> peptides = sequence.extractPeptides();
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
		U.p("Welcome to JavaGFS");
		U.p("Genome-based peptide fingerprint scanning");
		U.p("Developed 2010 by the Giddings Lab");
	}
	

}
