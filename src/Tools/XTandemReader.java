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
 * Jainab gave me some X!Tandem reports.  This reads them in and interprets the data
 * @author Brian Risk
 *
 */
public class XTandemReader {
	
	
	//the suffix that the files have
	static String suffix = ".xml";
	
	//what surrounds the spectrum file name

	static String spectrumStop = "'\">";
	static String eStart = "expect=\"";
	static String seqStart = "seq=\"";
	
	public static void main(String args[]) {
		
		ArrayList<TestSet> testSets = new ArrayList<TestSet>();
		
		testSets.add(getTestSet(
			"/Users/risk2/PeppyOverflow/reports - saved/X!Tandem reports/EColi",
			"/Users/risk2/PeppyOverflow/tests/",
			"ecoli",
			Color.red,
			"label=\"models from '/Users/khatun/PeppyToRunproteinDatabase/TestData/ecoli/spectra/"
		));
		
		testSets.add(getTestSet(
			"/Users/risk2/PeppyOverflow/reports - saved/X!Tandem reports/Kapp",
			"/Users/risk2/PeppyOverflow/tests/",
			"human",
			Color.blue,
			"label=\"models from '/Users/khatun/PeppyToRunproteinDatabase/TestData/kapp/spectra/"
		));
		
		testSets.add(getTestSet(
			"/Users/risk2/PeppyOverflow/reports - saved/X!Tandem reports/Aurum",
			"/Users/risk2/PeppyOverflow/tests/",
			"aurum",
			Color.green,
			"label=\"models from '/Users/khatun/PeppyToRunproteinDatabase/TestData/aurum/spectra/"
		));
		
		
		
		for (TestSet test: testSets) {
			test.cleanMatches();
			test.calculateStastics();
		}
		
		ValidationReport vr = new ValidationReport(testSets);
		vr.createReport();
		U.p("done");
	}
	
	private static TestSet getTestSet(String reportLocation, String testLocation, String testName, Color color, String spectrumStart) {
		
		/* get files from report folder */
		File reportFolder = new File(reportLocation);
		File [] reportFiles = reportFolder.listFiles();
		
		/* where we hold matches */
		ArrayList<Match> matches = new ArrayList<Match>();
		Match match;
		
		/* load the spectra from test set so we can get their spectrum IDs */
		ArrayList<Spectrum> spectra = Spectrum.loadSpectraFromFolder(testLocation + testName + "/spectra");
			
		for (int fileIndex = 0; fileIndex < reportFiles.length; fileIndex++) {
			if (!reportFiles[fileIndex].getName().toLowerCase().endsWith(suffix)) continue;
			match = extractHitsFromFile(reportFiles[fileIndex], spectrumStart, testLocation, testName);
			if (match != null) {
				
				/* set the proper spectrum ID */
				Spectrum matchSpectrum = match.getSpectrum();
				for (Spectrum spectrum: spectra) {
					if (spectrum.getFile().getName().equals(matchSpectrum.getFile().getName())) {
						matchSpectrum.setId(spectrum.getId());
					}
				}
				
				/* add the match */
				matches.add(match);
			}
		}
		U.p("found this many matches: " + matches.size());
		
		return new TestSet(testLocation, testName, matches, color);
		
	}
	
	private static Match_Blank extractHitsFromFile(File file, String spectrumStart, String testLocation, String testName) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			Spectrum spectrum = null;
			String line = br.readLine();
			String nugget;
			int subIndexStart, subIndexStop;
			while (line != null) {
				//the line which contains the spectrum
				if (line.startsWith("<bioml")) {
					subIndexStart = line.indexOf(spectrumStart) + spectrumStart.length();
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
