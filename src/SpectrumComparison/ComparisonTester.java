package SpectrumComparison;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;

import Peppy.Peptide;
import Peppy.Properties;
import Peppy.Spectrum;
import Utilities.U;

public class ComparisonTester {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new ComparisonTester("USP", 0.3);
		U.p("done");	
	}
	
	public ComparisonTester(String species, double delta) {
		
		Properties.highIntensityCleaning = true;
		
		ArrayList<SpectrumPeptidePair> spectrumPeptidePairs = new ArrayList<SpectrumPeptidePair>();

		//go through each file in our peptides folder
		File peptideFolder = new File("/Users/risk2/PeppyOverflow/tests/" + species + "/peptides");
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
			File spectrumFile = new File("/Users/risk2/PeppyOverflow/tests/" + species + "/spectra/" + peptideFiles[peptideFileIndex].getName());
			if (!spectrumFile.exists()) continue;
			
			/*
			 * Here we include the peptide we know to be correct and see how it would score	
			 */
			Peptide peptide = new Peptide(peptideString);
			//when we load the spectrum we set the flag to    normalize the peaks
			Spectrum spectrum = new Spectrum(spectrumFile, false);
			spectrumPeptidePairs.add(new SpectrumPeptidePair(spectrum, peptide));		
		}
		
		ArrayList<SpectrumComparison> spectrumComparisons = new ArrayList<SpectrumComparison>();
		for (int i = 0; i < spectrumPeptidePairs.size() - 1; i++) {
			SpectrumPeptidePair spp1 = spectrumPeptidePairs.get(i);
			for (int j = i + 1; j < spectrumPeptidePairs.size(); j++) {
				SpectrumPeptidePair spp2 = spectrumPeptidePairs.get(j);
				double precursorDifference = Math.abs(spp1.getSpectrum().getPrecursorMass() - spp2.getSpectrum().getPrecursorMass());
				if (precursorDifference < 2.0 ) {
					SpectrumComparison comparison = new SpectrumComparison(spp1,spp2, delta);
					comparison.computeDistance();
					spectrumComparisons.add(comparison);
				}
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
		
		//Generate a visual report, please!
		U.p("generating report");
		generateReport(spectrumComparisons);
		
	}
	
	public static void generateReport( ArrayList<SpectrumComparison> spectrumComparisons) {
		//intro!
		U.p("printing out the top matches");
		
		Collections.sort(spectrumComparisons);
		
		int width = 1000;
		int height = 200;
		
		File testFileFolder = new File(Properties.reportDirectory, "spectrumComparison");
		testFileFolder.mkdirs();
		File testFile = new File(testFileFolder, "index.html");
		File imagesFolder = new File(testFileFolder, "images");
		imagesFolder.mkdir();
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(testFile)));
			
			//print the header of the file
			pw.println("<html><body>");
			
			for (int i = 0; i < 50; i++) {
				SpectrumComparison comparison = spectrumComparisons.get(i);
				
				//creating the image of our match
				File ourMatchFile = new File (imagesFolder, i + "a.jpg");
				SpectralVisualizer.SpectralVisualizer.drawSpectrum(comparison.getSpectrumA(), width, height, ourMatchFile, false);
				
				//creating the image of their "true" match.  Whatever!
				File trueMatchFile = new File (imagesFolder, i + "b.jpg");
				SpectralVisualizer.SpectralVisualizer.drawSpectrum(comparison.getSpectrumB(), width, height, trueMatchFile, false);
				
				
				//print int out
				pw.print(comparison.getSpectrumA().getFile().getName());
				pw.print(";");
				pw.print(comparison.getSpectrumPeptidePairA().getPeptide().getAcidSequenceString());
				pw.print("<p>");
				pw.print(comparison.getSpectrumB().getFile().getName());
				pw.print(";");
				pw.print(comparison.getSpectrumPeptidePairB().getPeptide().getAcidSequenceString());
				pw.print("<p>");
				pw.print("<img src=\"");
				pw.print(imagesFolder.getName() + "/" + ourMatchFile.getName());
				pw.print("\"><br>");
				pw.print("<img src=\"");
				pw.print(imagesFolder.getName() + "/" + trueMatchFile.getName());
				pw.print("\"><br>");
				pw.print("<hr>");
				pw.print("<p>");
				pw.println();
				pw.println();
			}
		
			
			//print the footer and close out
			pw.println("</body></html>");
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}

}
