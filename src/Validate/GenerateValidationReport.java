package Validate;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;

import Peppy.JavaGFS;
import Peppy.Peptide;
import Peppy.ProteinDigestion;
import Peppy.Spectrum;
import Peppy.SpectrumPeptideMatch;
import Utilities.U;

public class GenerateValidationReport {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		File validationFolder = new File("Validation");
		validationFolder.mkdir();
		
		File indexFile = new File(validationFolder, "index.html");
		
		//load the swis prot protein database
		ArrayList<Peptide> peptides = ProteinDigestion.getPeptidesFromProteinFile(new File("tests/databases/uniprot_sprot.fasta"));
		
		//set up which tests we will perfrom
		ArrayList<String> tests = new ArrayList<String>();
		tests.add("human");
		tests.add("aurum");
		tests.add("ecoli");
		
		//get the matches for each of our tests
		for (int i = 0; i < tests.size(); i++) {
			String test = tests.get(i);
			U.p("Getting matches for: " + test);
			
			//load spectra for this test
			ArrayList<Spectrum> spectra = Spectrum.loadSpectraFromFolder("tests/" + test + "/spectra");
			U.p("loaded " +spectra.size() + " spectra.");
			
			//get the matches
			ArrayList<SpectrumPeptideMatch> matches = JavaGFS.asynchronousDigestion(peptides, spectra, null);
			
			//See which of these matches are true
			ArrayList<MatchContainer> testedMatches = new ArrayList<MatchContainer>(matches.size());
			for (SpectrumPeptideMatch match: matches) {
				testedMatches.add(new MatchContainer(match));
			}
			
			//Sort by e value
			Collections.sort(testedMatches);
		}
		 
		//get all the matches
			
		

	}
	
	public static void generatePrecisionRecallCurve(ArrayList<SpectrumPeptideMatch> matches) {
		
	}

}
