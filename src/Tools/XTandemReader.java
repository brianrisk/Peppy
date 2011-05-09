package Tools;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import Peppy.Match;
import Peppy.Match_Blank;
import Peppy.Peptide;
import Peppy.Spectrum;
import Utilities.U;
import Validate.TestSet;
import Validate.ValidationReport;

/**
 * Jainab gave me some X!Tandem reports.  This reads them in and interprets the data
 * @author Brian Risk
 *
 */
public class XTandemReader {
	
	//Where our report files are located
	static String reportLocation = "reports/X!Tandem reports/EColi";
//	static String reportLocation = "reports/X!Tandem reports/Kapp";
//	static String reportLocation = "reports/X!Tandem reports/Aurum";
	
	//the suffix that the files have
	static String suffix = ".xml";
	
	//test name (from ValidationReport)
	static String testName = "ecoli";
//	static String testName = "human";
//	static String testName = "aurum";
	
	static Color color = Color.red;
//	static Color color = Color.green;
//	static Color color = Color.blue;
	
	//test spectra location
	static String testLocation = "/Users/risk2/PeppyOverflow/tests/";
	
	//what surrounds the spectrum file name
	static String specrumStart = "label=\"models from '/Users/khatun/PeppyToRunproteinDatabase/TestData/ecoli/spectra/";
//	static String specrumStart = "label=\"models from '/Users/khatun/PeppyToRunproteinDatabase/TestData/kapp/spectra/";
//	static String specrumStart = "label=\"models from '/Users/khatun/PeppyToRunproteinDatabase/TestData/aurum/spectra/";
	static String spectrumStop = "'\">";
	static String eStart = "expect=\"";
	static String seqStart = "seq=\"";
	
	public static void main(String args[]) {
		File reportFolder = new File(reportLocation);
		File [] reportFiles = reportFolder.listFiles();
		ArrayList<Match> matches = new ArrayList<Match>();
		Match match;
		for (int fileIndex = 0; fileIndex < reportFiles.length; fileIndex++) {
			if (!reportFiles[fileIndex].getName().toLowerCase().endsWith(suffix)) continue;
			match = extractHitsFromFile(reportFiles[fileIndex]);
			if (match != null) matches.add(match);
		}
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
	
	private static Match_Blank extractHitsFromFile(File file) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			Spectrum spectrum = null;
			String line = br.readLine();
			String nugget;
			int subIndexStart, subIndexStop;
			while (line != null) {
				//the line which contains the spectrum
				if (line.startsWith("<bioml")) {
					subIndexStart = line.indexOf(specrumStart) + specrumStart.length();
					subIndexStop = line.indexOf(spectrumStop, subIndexStart);
					nugget = line.substring(subIndexStart, subIndexStop);
					spectrum = new Spectrum(testLocation + testName + "/spectra/" + nugget);

				}
				
				//line which contains match data like peptide, evalue, score
				if (line.startsWith("<domain")) {
					//get the evalue
					subIndexStart = line.indexOf(eStart) + eStart.length();
					subIndexStop = line.indexOf("\"", subIndexStart);
					nugget = line.substring(subIndexStart, subIndexStop);
					double evalue = Double.parseDouble(nugget);
					
					//get the acid sequence
					subIndexStart = line.indexOf(seqStart) + seqStart.length();
					subIndexStop = line.indexOf("\"", subIndexStart);
					nugget = line.substring(subIndexStart, subIndexStop);
					
					return new Match_Blank(spectrum, new Peptide(nugget), 0.0, evalue);
					
					
				}
				
				line=br.readLine();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}

}
