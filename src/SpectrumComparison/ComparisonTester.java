package SpectrumComparison;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

import Peppy.Peptide;
import Peppy.Spectrum;
import Utilities.U;

public class ComparisonTester {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new ComparisonTester("human", 0.3);
		U.p("done");	
	}
	
	public ComparisonTester(String species, double delta) {
		ArrayList<SpectrumPeptidePair> spectrumPeptidePairs = new ArrayList<SpectrumPeptidePair>();
		//Load our spectra
		ArrayList<Spectrum> spectra = Spectrum.loadSpectraFromFolder("tests/" + species + "/spectra");
		//go through each file in our peptides folder
		File peptideFolder = new File("tests/" + species + "/peptides");
		File [] peptideFiles = peptideFolder.listFiles();
		for (int peptideFileIndex = 0; peptideFileIndex < peptideFiles.length; peptideFileIndex++) {
			//only want visible files
			if (peptideFiles[peptideFileIndex].isHidden()) continue;
			
			//only want valid file types
			String fileName = peptideFiles[peptideFileIndex].getName();
			if (!fileName.endsWith(".dta") && !fileName.endsWith(".pkl") && !fileName.endsWith(".txt")) continue;
			
			//load in our string line
			String peptideString = "";
			try {
				BufferedReader br = new BufferedReader(new FileReader(peptideFiles[peptideFileIndex]));
				//read the first line;
				peptideString = br.readLine();
				//close;
				br.close();
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
			//Okay, we've got a valid peptide file
			if (peptideString == null) continue;
			peptideString = peptideString.trim();
			if (peptideString.equals("")) continue;
			
			//make sure this peptide has a corresponding spectrum file
			File spectrumFile = new File("tests/" + species + "/spectra/" + peptideFiles[peptideFileIndex].getName());
			if (!spectrumFile.exists()) continue;
			
			/*
			 * Here we include the peptide we know to be correct and see how it would score	
			 */
			Peptide peptide = new Peptide(peptideString);
			//This assumes only one spectrum per file, but that should be the case with these test cases.
			Spectrum spectrum = Spectrum.loadSpectra(spectrumFile).get(0);
			
			spectrumPeptidePairs.add(new SpectrumPeptidePair(spectrum, peptide));		
		}
		
		ArrayList<SpectrumComparison> spectrumComparisons = new ArrayList<SpectrumComparison>();
		for (int i = 0; i < spectrumPeptidePairs.size() - 1; i++) {
//			U.p(i + 1 + " of " + spectrumPeptidePairs.size());
			for (int j = i + 1; j < spectrumPeptidePairs.size(); j++) {
//				U.p("sub: " + j);
				SpectrumComparison comparison = new SpectrumComparison(spectrumPeptidePairs.get(i), spectrumPeptidePairs.get(j), delta);
				comparison.computeDistance();
				spectrumComparisons.add(comparison);
			}
		}
		U.p("There are this many compariosns: " + spectrumComparisons.size());
		Collections.sort(spectrumComparisons);
		for (int i = 0; i < spectrumComparisons.size(); i++) {
			if (!spectrumComparisons.get(i).isEqual()) {
				U.p("the cutoff distance: " + spectrumComparisons.get(i).getDistance());
				U.p("We found this many matches: " + i);
				break;
			}
		}
		
		for (int i = 0; i < 10; i ++) {
			U.p(spectrumComparisons.get(i));
		}
		
	}

}
