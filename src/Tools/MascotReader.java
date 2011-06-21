package Tools;


import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

import Peppy.Match;
import Peppy.Match_Blank;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.Spectrum;
import Utilities.U;
import Validate.TestSet;
import Validate.ValidationReport;

/**
 * This interprets the files Yanbao gave me
 * @author Brian Risk
 *
 */
public class MascotReader {
	
	
	public static void main(String args[]) {
		Properties.maximumNumberOfMatchesForASpectrum = 1;
		
		ArrayList<TestSet> testSets = new ArrayList<TestSet>();
		testSets.add(getTestSet(
				"reports/Mascot/ecoli.txt",
				"/Users/risk2/PeppyOverflow/tests/",
				"ecoli",
				Color.red
		));
		
		testSets.add(getTestSet(
				"reports/Mascot/kapp.txt",
				"/Users/risk2/PeppyOverflow/tests/",
				"human",
				Color.blue
		));
		
		testSets.add(getTestSet(
				"reports/Mascot/aurum.txt",
				"/Users/risk2/PeppyOverflow/tests/",
				"aurum",
				Color.green
		));
		
		for (TestSet test: testSets) {
			test.calculateStastics();
		}
		
		ValidationReport vr = new ValidationReport(testSets);
		vr.createReport();
		U.p("done");
	}
	
	private static TestSet getTestSet(String reportLocation, String testLocation, String testName, Color color) {
		//where our matches are
		File matchesFile = new File(reportLocation);
		
		//the spectra of our tests
		ArrayList<Spectrum> spectra = Spectrum.loadSpectraFromFolder(new File(testLocation + testName + "/" + "spectra/"));
		
		//set spectra id according to mass
		Collections.sort(spectra);
		for (int i = 0; i < spectra.size(); i++) {
			spectra.get(i).setId(i);
		}
		
		//get the matches that correspond with the spectra
		ArrayList<Match> matches = extractMatchesFromFile(matchesFile, spectra);
		
		return new TestSet(testLocation, testName, matches, color);
	}
	
	private static ArrayList<Match> extractMatchesFromFile(File file, ArrayList<Spectrum> spectra) {
		ArrayList<Match> matches = new ArrayList<Match>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			String line = br.readLine();
			
			int id;
			double eValue;
			String peptide;
			while (line != null) {
				String [] chunks = line.split("\t");
				id = Integer.parseInt(chunks[0]) - 1;
				eValue = Double.parseDouble(chunks[1]);
				peptide = chunks[2];
				
				Match_Blank match = new Match_Blank(spectra.get(id), new Peptide(peptide), 0.0, eValue);
				
				//look for other match with this id
				//delete if one exists with higher e value
				boolean addMatch = true;
				for (int i = 0; i < matches.size(); i++) {
					if (matches.get(i).getSpectrum().getId() == id) {
						if (matches.get(i).getEValue() < eValue) {
							addMatch = false;
						} else {
							if (matches.get(i).getEValue() > eValue) {
								matches.remove(i);
							}
						}
					}
				}
				
				if (addMatch) matches.add(match);

				line=br.readLine();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return matches;
	}




}
