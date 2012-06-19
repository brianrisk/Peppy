package Tools;


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
import Peppy.SpectrumLoader;
import Peppy.U;
import Validate.TestSet;
import Validate.ValidationReport;

/**
 * This interprets the files Yanbao gave me
 * @author Brian Risk
 *
 */
public class MascotReader {
	
	
	public static void main(String args[]) {
		
		ArrayList<TestSet> testSets = new ArrayList<TestSet>();
		testSets.add(getTestSet(
				"/Users/risk2/PeppyData/reports - saved/Mascot/ecoli.txt",
				"/Users/risk2/PeppyData/tests/",
				"ecoli"
		));
		
		testSets.add(getTestSet(
				"/Users/risk2/PeppyData/reports - saved/Mascot/kapp.txt",
				"/Users/risk2/PeppyData/tests/",
				"human"
		));
		
		testSets.add(getTestSet(
				"/Users/risk2/PeppyData/reports - saved/Mascot/aurum.txt",
				"/Users/risk2/PeppyData/tests/",
				"aurum"
		));
		
		for (TestSet test: testSets) {
			test.calculateStastics();
		}
		
		ValidationReport vr = new ValidationReport(testSets);
		Properties.scoringMethodName = "Mascot";
		vr.createReport();
		vr.createResultsFiles();
		U.p("done");
	}
	
	private static TestSet getTestSet(String reportLocation, String testLocation, String testName) {
		//where our matches are
		File matchesFile = new File(reportLocation);
		
		//the spectra of our tests
		ArrayList<Spectrum> spectra = SpectrumLoader.loadSpectraFromFolder(new File(testLocation + testName + "/" + "spectra/"));
		
		//set spectra id according to mass
		Collections.sort(spectra);
		for (int i = 0; i < spectra.size(); i++) {
			spectra.get(i).setId(i);
		}
		
		//get the matches that correspond with the spectra
		ArrayList<Match> matches = extractMatchesFromFile(matchesFile, spectra);
		
		return new TestSet(testLocation, testName, matches, spectra);
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
				
				Match_Blank match = new Match_Blank(spectra.get(id), new Peptide(peptide), eValue);
				matches.add(match);

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
