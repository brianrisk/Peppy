package Validate;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import Peppy.Peptide;
import Peppy.Properties;
import Peppy.Sequence_Protein;
import Peppy.Spectrum;
import Utilities.U;

/**
 * Test sets often contain matches that will not be found with
 * normal search parameters.  This class takes a specified
 * test set and produces a subset of the test which has vetted matches
 * @author Brian Risk
 *
 */
public class TestSetReduce {
	
	private static final String testDirectoryName = "/Users/risk2/PeppyOverflow/tests/";
	
	public static void main(String args[]) {
		
		//load our peptide database
		U.p("loading database");
		Properties.numberOfMissedCleavages = 2;
		File databaseFile = new File("/Users/risk2/PeppyOverflow/tests/databases/uniprot_sprot.fasta");
		Sequence_Protein sequence = new Sequence_Protein(databaseFile);
		ArrayList<Peptide> peptides = sequence.extractAllPeptides(false);
		
		//specify test set
		String testName;
//		String [] testNames = {"ecoli", "human", "aurum", "USP"};
		String [] testNames = {"human"};
		for (int testIndex = 0; testIndex < testNames.length; testIndex++) {
			testName = testNames[testIndex];
			U.p("doing: " + testName);
		
		
			//create folders where we will store our output
			File newTestFolder = new File (testDirectoryName, testName + "-vetted");
			newTestFolder.mkdirs();
			File newPeptidesFolder = new File(newTestFolder, "peptides");
			newPeptidesFolder.mkdirs();
			File newSpectraFolder = new File(newTestFolder, "spectra");
			newSpectraFolder.mkdirs();
			
			
			
			//load spectra
			ArrayList<Spectrum> spectra = Spectrum.loadSpectraFromFolder(testDirectoryName + testName + "/spectra");
			U.p("original spectra size: " + spectra.size());
			
			//eliminate duplicate spectra
			int iMax = spectra.size() - 1;
			int jMax = spectra.size();
			Spectrum spectrumI, spectrumJ;
			for (int i = 0; i < iMax; i++) {
				spectrumI = spectra.get(i);
				for (int j = i + 1; j < jMax; j++) {
					spectrumJ = spectra.get(j);
					if (spectrumI.getMD5().equals(spectrumJ.getMD5())) {
						spectra.remove(j);
						j--;
						jMax--;
						iMax--;
					}
				}
			}
			U.p("spectra size without duplicates: " + spectra.size());
			
			//loop through spectra and add those with valid peptide files
			for (Spectrum spectrum: spectra) {
				//find the file for the correct peptide
				File spectrumFile = spectrum.getFile();
				File testFolder = spectrumFile.getParentFile().getParentFile();
				File peptideFolder = new File(testFolder, "peptides");
				File peptideFile = new File(peptideFolder, spectrumFile.getName());
				
				//load in the correct peptide string
				String correctAcidSequence = "";
				try {
					BufferedReader br = new BufferedReader(new FileReader(peptideFile));
					//read the first line;
					correctAcidSequence = br.readLine();
					//close;
					br.close();
				} catch (FileNotFoundException e) {
					continue;
				} catch (IOException e) {
					continue;
				}
				
				//testing that we've got a valid peptide file
				if (correctAcidSequence == null) continue;
				correctAcidSequence = correctAcidSequence.trim();
				if (correctAcidSequence.equals("")) continue;
				
				if (testName.equals("ecoli")) {
					if (correctAcidSequence.equals("KGEMNFDVVIASPDAMR")) continue;
					if (correctAcidSequence.equals("EQIIFPEIDYDKVDR")) continue;
				}
				
				//create a new peptide from out loaded string
				Peptide peptide = new Peptide(correctAcidSequence);
				
				//see that the peptide matches the minimum weight
				if (peptide.getMass() < Properties.peptideMassThreshold) continue;
				
				//see that the difference between the predicted mass and the precursor mass is within tolerance
				if (Math.abs(peptide.getMass() - spectrum.getMass()) > Properties.precursorTolerance) {
//					U.p("Peptided outside of mass error range: " + peptide.getAcidSequenceString());
					continue;
				}
				
				//see that the peptide is in the database
				if (ValidationReport.isPeptidePresentInList( peptide, peptides) < 0) {
//					U.p("Peptided not in database: " + peptide.getAcidSequenceString());
					continue;
				}
				
				//if all these steps have passed, it must be good.  save files
				File newSpectrumFile = new File(newSpectraFolder, spectrumFile.getName());
				File newPeptideFile = new File(newPeptidesFolder, peptideFile.getName());
				U.copyfile(spectrumFile, newSpectrumFile);
				U.copyfile(peptideFile, newPeptideFile);
				
			}
		}
		
		U.p("done");
		
	}

}
