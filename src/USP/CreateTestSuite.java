package USP;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import Peppy.Match;
import Peppy.MatchConstructor;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.ScoringThreadServer;
import Peppy.Sequence_Protein;
import Peppy.Spectrum;
import Utilities.U;

public class CreateTestSuite {

	/**
	 * This uses the USP data to create a test suite similar to that of the
	 * Aurum data set.  This method is not perfect, but it may be good enough
	 * for initial testing and optimization.
	 * 
	 * The process is this:  We sue the small database of peptides gleaned from
	 * ExtractListFromDatabase and then runs the spectral set against it using
	 * MSMSFit.  The idea is that since the database is comparatively small that
	 * the top match will be the correct match as the odds of a false positive
	 * in such a small database are slim.  The peptide sequences are then stored
	 * in a way similar to the way the Aurum set -- each peptide sequence has
	 * its own text file which is named the same as the spectrum it is paired
	 * with.
	 * @param args
	 */
	public static void main(String[] args) {
		//Opening message
		U.p("Creating the necessary files for our USP test suite");
//		createSuiteFromUSPDatabase();
		createSuiteFromTopX(10);
//		createSuiteWithClosePrecursors();
		
		U.p("Done!");
	}
	
	/**
	 * finds X peptide matches for each spectrum.  Goes through and finds which
	 * of these matches have peptides from the USP 50 database.  Those go into our
	 * test set.
	 */
	public static void createSuiteFromTopX(int numberToKeepPerSpectrum) {
		//load the peptides from the database
		Properties.numberOfMissedCleavages = 1;
		//load the correct peptide set
		Sequence_Protein sequenceCorrect = new Sequence_Protein(new File("/Users/risk2/PeppyOverflow/USP/extracted-proteins.txt"));
		ArrayList<Peptide> correctPeptides = sequenceCorrect.extractAllPeptides(false);
		
		//load full database
		Sequence_Protein sequence = new Sequence_Protein(new File("/Users/risk2/PeppyOverflow/tests/databases/uniprot_sprot.fasta"));
		ArrayList<Peptide> peptides = sequence.extractAllPeptides(false);
		
		//Load our spectra
		Properties.spectraDirectoryOrFile = new File("/Users/risk2/PeppyOverflow/spectra USP");
		ArrayList<Spectrum> spectra = Spectrum.loadSpectra();
		U.p("loaded " +spectra.size() + " spectra.");

		//restrict the precursors since they should be pretty accurate
		Properties.precursorTolerance = 0.7;
		
		//Get the  matches
		Properties.maximumNumberOfMatchesForASpectrum = numberToKeepPerSpectrum;
		Properties.scoringMethodName = "Peppy.Match_IMP";
		Properties.matchConstructor = new MatchConstructor(Properties.scoringMethodName);
		ArrayList<Match> matches  = (new ScoringThreadServer(peptides, spectra)).getMatches();
		
		//save to appropriate files
		File peptideDir = new File("/Users/risk2/PeppyOverflow/tests/USP top " + numberToKeepPerSpectrum + "/peptides/");
		peptideDir.mkdirs();
		File spectraFolder =  new File("/Users/risk2/PeppyOverflow/tests/USP top " + numberToKeepPerSpectrum + "/spectra/");
		spectraFolder.mkdirs();
		for (Match match: matches) {
			try {
				for (Peptide peptide: correctPeptides) {
					if (match.getPeptide().equals(peptide)) {
						File peptideFile = new File(peptideDir, match.getSpectrum().getFile().getName());
						PrintWriter pw;
						pw = new PrintWriter(new BufferedWriter(new FileWriter(peptideFile)));
						pw.println(match.getPeptide().getAcidSequenceString());
						pw.flush();
						pw.close();
						
						//copy the spectrum file
						File spectraFile = new File(spectraFolder, match.getSpectrum().getFile().getName());
						spectraFolder.mkdirs();
						U.copyfile(match.getSpectrum().getFile(), spectraFile);
					}
				}					
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	public static void createSuiteWithClosePrecursors() {
		//set our properties
		Properties.spectraDirectoryOrFile = new File("/Users/risk2/PeppyOverflow/spectra USP");
		Properties.maximumNumberOfMatchesForASpectrum = 1;
		
		//load the peptides from the database
		ArrayList<Peptide> peptides = new ArrayList<Peptide>();
		File uniqePeptidesFile = new File("/Users/risk2/PeppyOverflow/USP/peptide-unique-mass-list.txt");
		try {
			BufferedReader br = new BufferedReader(new FileReader(uniqePeptidesFile));
			String peptideLine;
			peptideLine = br.readLine();
			while (peptideLine != null) {
				peptideLine = peptideLine.trim();
				if (!peptideLine.equals("")) {
					peptides.add(new Peptide(peptideLine));
				}
				peptideLine = br.readLine();
			}
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		U.p("laded " + peptides.size() + " peptides");
		
		//Load our spectra
		ArrayList<Spectrum> spectra = Spectrum.loadSpectra();
		U.p("loaded " +spectra.size() + " spectra.");

		
		//make test directory
		File spectraFolder =  new File("/Users/risk2/PeppyOverflow/tests/USP-0.06 precursor tolerance/spectra/");
		File peptideFolder =  new File("/Users/risk2/PeppyOverflow/tests/USP-0.06 precursor tolerance/peptides/");
		spectraFolder.mkdirs();
		peptideFolder.mkdirs();
	
		
		//go through each spectrum, find those which have close precursor to peptide
		for (Spectrum spectrum: spectra) {
			
			for (Peptide peptide: peptides) {
				
				if (Math.abs(spectrum.getMass() - peptide.getMass()) < 0.06) {
					//save the peptide file
					File peptideFile = new File(peptideFolder, spectrum.getFile().getName());
					PrintWriter pw;
					try {
						pw = new PrintWriter(new BufferedWriter(new FileWriter(peptideFile)));
						pw.println(peptide.getAcidSequenceString());
						pw.flush();
						pw.close();
					} catch (IOException e) {
						e.printStackTrace();
					}
					
					//copy the spectrum file
					File spectraFile = new File(spectraFolder, spectrum.getFile().getName());
					U.copyfile(spectrum.getFile(), spectraFile);
					break;
				}
			}
		}
		
	}

}
