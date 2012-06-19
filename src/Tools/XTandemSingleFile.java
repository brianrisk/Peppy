package Tools;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

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
 * For analyzing the files John created with X!Tandem
 * @author Brian Risk
 *
 */
public class XTandemSingleFile {
	
	
	public static void main(String args[]) {
		Properties.scoringMethodName = "XTandem";
		
		ArrayList<TestSet> testSets = new ArrayList<TestSet>();
		
		testSets.add(getTestSet(
			"/Users/risk2/PeppyData/reports - saved/X!Tandem John/output_USP-top-10.2012_05_09_14_46_03.t.xml",
			"/Users/risk2/PeppyData/tests/",
			"USP top 10"
		));
		
		testSets.add(getTestSet(
			"/Users/risk2/PeppyData/reports - saved/X!Tandem John/output_human.2012_05_09_13_25_02.t.xml",
			"/Users/risk2/PeppyData/tests/",
			"human"
		));
		
		testSets.add(getTestSet(
			"/Users/risk2/PeppyData/reports - saved/X!Tandem John/output_aurum.2012_05_09_14_36_22.t.xml",
			"/Users/risk2/PeppyData/tests/",
			"aurum"
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
		
		/* get files from report folder */
		File reportFile = new File(reportFileName);
		
		/* where we hold matches */
		ArrayList<Match> matches = new ArrayList<Match>();
		Match match;
		
		/* load the spectra from test set so we can get their spectrum IDs */
		ArrayList<Spectrum> spectra = SpectrumLoader.loadSpectraFromFolder(testLocation + testName + "/spectra");
		
		/* load and parse the XML */
		try {
			DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
			DocumentBuilder documentBuilder = factory.newDocumentBuilder();
			Document doc = documentBuilder.parse (reportFile);

            // normalize text representation
            doc.getDocumentElement ().normalize ();
            System.out.println ("Root element of the doc is " + doc.getDocumentElement().getNodeName());

            NodeList listOfNodes = doc.getElementsByTagName("domain");
            
            String spectrum;
            String eValue;
            String peptide;
            int spectrumNumber;
            for (int i = 0; i < listOfNodes.getLength(); i++) {
            	Node node = listOfNodes.item(i);
            	spectrum = ((Element)node).getAttribute("id");
            	eValue  = ((Element)node).getAttribute("expect");
            	peptide  = ((Element)node).getAttribute("seq");
            	
            	spectrumNumber = Integer.parseInt(spectrum.split("\\.")[0]) - 1;
            	
            	/* do we have to take the neg log of the e value here? */
            	match = new Match_Blank(spectra.get(spectrumNumber), new Peptide(peptide), Double.parseDouble(eValue));
            	matches.add(match);
            }

		} catch (ParserConfigurationException e) {
			e.printStackTrace();
		} catch (SAXException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
			
		
		
		U.p("found this many matches: " + matches.size());
		
		return new TestSet(testLocation, testName, matches, spectra);
		
	}
	
	


}
