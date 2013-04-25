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
 * built do interpret the second batch of mascot files that Yanbao created
 * @author Brian Risk
 *
 */
public class MascotReader2 {
	
	public static void main2x(String args[]) {
		File reportsDir = new File ("/Users/risk2/PeppyData/reports - saved/Mascot 2/");
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(reportsDir, "aurum-unassigned.csv")));
			String line = br.readLine();
			while (line!=null) {
				String [] chunks = line.split(",");
				if (chunks.length < 25) {
					U.p(line);
//					U.p(chunks[12] + "\t" + chunks[24] + "\t" + chunks[27]);
				}
				line = br.readLine();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String args[]) {
		Properties.scoringMethodName = "Mascot";
		ArrayList<TestSet> testSets = new ArrayList<TestSet>();
		
		testSets.add(getTestSet(
				"/Users/risk2/PeppyData/reports - saved/mascot-reports/kapp-unassigned2.csv",
				"/Users/risk2/PeppyData/tests/",
				"human"
		));
		
		testSets.add(getTestSet(
				"/Users/risk2/PeppyData/reports - saved/mascot-reports/aurum-unassigned2.csv",
				"/Users/risk2/PeppyData/tests/",
				"aurum"
		));
		
		testSets.add(getTestSet(
				"/Users/risk2/PeppyData/reports - saved/mascot-reports/ups-unassigned2.csv",
				"/Users/risk2/PeppyData/tests/",
				"USP top 10"
		));
		
		for (TestSet test: testSets) {
			test.calculateStastics();
		}
		
		ValidationReport vr = new ValidationReport(testSets);
		vr.createReport();
		vr.createResultsFiles();
		U.p("done");
	}
	
	private static TestSet getTestSet(String reportFileName, String testLocation, String testName) {
		
		/* where we hold matches */
		ArrayList<Match> matches = new ArrayList<Match>();
		Match match;
		
		/* load the spectra from test set so we can get their spectrum IDs */
		ArrayList<Spectrum> spectra = SpectrumLoader.loadSpectra(testLocation + testName + "/spectra");
		
		/* sorting by mass as that's how they are sorted in the mascot results (thus giving them new IDs...) */
		Collections.sort(spectra);
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(reportFileName)));
			String line = br.readLine();
			line = br.readLine();
			while (line!=null) {
				String [] chunks = line.split(",");
				if (chunks.length < 3) {
					line = br.readLine();
					continue;
				}
				int spectrumNumber = Integer.parseInt(chunks[0]) - 1;
//				double score = Double.parseDouble(chunks[1]);
				double score = Double.parseDouble(chunks[3]);
				score *= -1;
				Peptide peptide = new Peptide(chunks[2]);
				match = new Match_Blank(spectra.get(spectrumNumber),peptide, score);
		    	matches.add(match);
				line = br.readLine();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} 

		U.p("found this many matches: " + matches.size());
		
		return new TestSet(testLocation, testName, matches, spectra);
		
	}

}
