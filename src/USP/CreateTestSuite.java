package USP;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import Peppy.Peppy;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.ProteinDigestion;
import Peppy.ScoringThreadServer;
import Peppy.Spectrum;
import Peppy.Match;
import Utilities.U;

public class CreateTestSuite {

	/**
	 * This uses the USP data to create a test suite similar to that of the
	 * Kapp data set.  This method is not perfect, but it may be good enough
	 * for initial testing and optimization.
	 * 
	 * The process is this:  We sue the small database of peptides gleaned from
	 * ExtractListFromDatabase and then runs the spectral set against it using
	 * MSMSFit.  The idea is that since the database is comparatively small that
	 * the top match will be the correct match as the odds of a false positive
	 * in such a small database are slim.  The peptide sequences are then stored
	 * in a way similar to the way the Kapp set -- each peptide sequence has
	 * its own text file which is named the same as the spectrum it is paired
	 * with.
	 * @param args
	 */
	public static void main(String[] args) {
		//Opening message
		U.p("Creating the necessary files for our USP test suite");
//		createSuiteFromUSPDatabase();
		createSuiteFromTopX(10);
		
		U.p("Done!");
	}
	
	/**
	 * finds X peptide matches for each spectrum.  Goes through and finds which
	 * of these matches have peptides from the USP 50 database.  Those go into our
	 * test set.
	 */
	public static void createSuiteFromTopX(int x) {
		//set our properties
		Properties.spectraDirectoryOrFile = new File("/Users/risk2/PeppyOverflow/spectra USP");
		Properties.maximumNumberOfMatchesForASpectrum = x;
		
		//Load our spectra
		ArrayList<Spectrum> spectra = Spectrum.loadSpectra();
		U.p("loaded " +spectra.size() + " spectra.");

		//load the peptides from the database
		ArrayList<Peptide> peptides = ProteinDigestion.getPeptidesFromFASTA(new File("/Users/risk2/PeppyOverflow/tests/databases/uniprot_sprot.fasta"));
		
		//load the correct peptide set
		ArrayList<Peptide> correctPeptides = ProteinDigestion.getPeptidesFromFASTA(new File("/Users/risk2/PeppyOverflow/USP/extracted-proteins.txt"));
		
		//restrict the precursors since they should be pretty accurate
		Properties.spectrumToPeptideMassError = 0.01;
		
		//Get the  matches
		ArrayList<Match> matches = (new ScoringThreadServer(peptides, spectra, null)).getMatches();
		
		//save to appropriate files
		File peptideDir = new File("/Users/risk2/PeppyOverflow/tests/USP/peptides/");
		peptideDir.mkdirs();
		File spectraFolder =  new File("/Users/risk2/PeppyOverflow/tests/USP/spectra/");
		spectraFolder.mkdirs();
		for (Match match: matches) {
			try {
				for (Peptide peptide: correctPeptides) {
					if (match.getPeptide().getAcidSequence().equals(peptide.getAcidSequence())) {
						File peptideFile = new File(peptideDir, match.getSpectrum().getFile().getName());
						PrintWriter pw;
						pw = new PrintWriter(new BufferedWriter(new FileWriter(peptideFile)));
						pw.println(match.getPeptide().getAcidSequence());
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
	
	public static void createSuiteFromUSPDatabase() {
		//set our properties
		Properties.spectraDirectoryOrFile = new File("/Users/risk2/PeppyOverflow/spectra USP");
		Properties.maximumNumberOfMatchesForASpectrum = 1;
		
		//Load our spectra
		ArrayList<Spectrum> spectra = Spectrum.loadSpectra();
		U.p("loaded " +spectra.size() + " spectra.");

		//load the peptides from the database
		ArrayList<Peptide> peptides = ProteinDigestion.getPeptidesFromFASTA(new File("/Users/risk2/PeppyOverflow/USP/extracted-proteins.txt"));
		
		//Get the matches
		ArrayList<Match> matches = (new ScoringThreadServer(peptides, spectra, null)).getMatches();
		
		//save to appropriate files
		for (Match match: matches) {
			try {
				File peptideFile = new File("/Users/risk2/PeppyOverflow/tests/USP/peptides/" + match.getSpectrum().getFile().getName());
				PrintWriter pw;
				pw = new PrintWriter(new BufferedWriter(new FileWriter(peptideFile)));
				pw.println(match.getPeptide().getAcidSequence());
				pw.flush();
				pw.close();
				
				//copy the spectrum file
				File spectraFolder =  new File("/Users/risk2/PeppyOverflow/tests/USP/spectra/");
				File spectraFile = new File(spectraFolder, match.getSpectrum().getFile().getName());
				spectraFolder.mkdirs();
				U.copyfile(match.getSpectrum().getFile(), spectraFile);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

}
