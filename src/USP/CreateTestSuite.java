package USP;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import Peppy.JavaGFS;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.ProteinDigestion;
import Peppy.Spectrum;
import Peppy.SpectrumPeptideMatch;
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
		
		//set our properties
		Properties.spectraDirectoryOrFile = new File("spectra USP");
		Properties.maximumNumberOfMatchesForASpectrum = 1;
		
		//Load our spectra
		ArrayList<Spectrum> spectra = Spectrum.loadSpectra();
		U.p("loaded " +spectra.size() + " spectra.");

		//load the peptides from the database
		ArrayList<Peptide> peptides = ProteinDigestion.getPeptidesFromProteinFile(new File("USP/extracted-proteins.txt"));
		
		//Get the matches
		ArrayList<SpectrumPeptideMatch> matches = JavaGFS.asynchronousDigestion(peptides, spectra, null);
		
		//save to appropriate files
		for (SpectrumPeptideMatch match: matches) {
			File peptideFile = new File("tests/USP/peptides/" + match.getSpectrum().getFile().getName());
			PrintWriter pw;
			try {
				pw = new PrintWriter(new BufferedWriter(new FileWriter(peptideFile)));
				pw.println(match.getPeptide().getAcidSequence());
				pw.flush();
				pw.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		U.p("Done!");

	}

}
