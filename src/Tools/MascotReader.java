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
	

	//Where our report files are located
//	static String reportLocation = "reports/Mascot/ecoli.txt";
//	static String reportLocation = "reports/Mascot/kapp.txt";
	static String reportLocation = "reports/Mascot/aurum.txt";
	
	//test name (from ValidationReport)
//	static String testName = "ecoli";
//	static String testName = "human";
	static String testName = "aurum";
	
//	static Color color = Color.red;
//	static Color color = Color.green;
	static Color color = Color.blue;
	
	//test spectra location
	static String testLocation = "/Users/risk2/PeppyOverflow/tests/";
	
	
	public static void main(String args[]) {

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
		
		U.p("found this many matches: " + matches.size());
		
		ArrayList<TestSet> testSets = new ArrayList<TestSet>();
		testSets.add(new TestSet(testLocation, testName, matches, color));
		
		for (TestSet test: testSets) {
			test.calculateStastics();
		}
		
		ValidationReport vr = new ValidationReport(testSets);
		vr.createReport();
		U.p("done");
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
							matches.remove(i);
							break;
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
