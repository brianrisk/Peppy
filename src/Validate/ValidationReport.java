	package Validate;

import java.awt.image.BufferedImage;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;

import javax.imageio.ImageIO;

import Math.MathFunctions;
import Peppy.Match;
import Peppy.MatchConstructor;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.Sequence;
import Peppy.Sequence_DNA;
import Peppy.Sequence_Protein;
import Peppy.U;

public class ValidationReport {
	
	public static ArrayList<TestSet> tests;
	public static File databaseFile = new File("");
	public static File reportFolder;
	public static PrintWriter indexWriter;

	public static int forwardsDatabaseSize = 0;
	public static int reverseDatabaseSize = 0;

	/**
	 * @param args
	 */
	public static void main(String[] args) {	
		Peppy.Peppy.init(args);
		setUp();
		addTests();
		U.startStopwatch();
//		createListOfPeptidesNotFoundInTheDatabase();
		searchTestSets();
		createReport();
		createResultsFiles();
//		createWrongReport();
		U.stopStopwatch();
		U.p("done.");
		
	}
	
	/**
	 * If a test set has been previously generated and loaded
	 * @param testSet
	 */
	public ValidationReport(ArrayList<TestSet> testSets) {
		tests = testSets;
		setUp();
	}
	
	
	private static void setUp() {
		System.setProperty("java.awt.headless", "true"); 
		
		/* Hello, world! */
		U.p("Are you ready for the food ball?  I mean: football.  I mean:  validation report");
		
		/* set up where we will save these reports */
		reportFolder =  new File(Properties.validationDirectory, "" + System.currentTimeMillis() + "/");
		reportFolder.mkdirs();
		File indexFile = new File(reportFolder, "index.html");
		try {
			indexWriter = new PrintWriter(new BufferedWriter(new FileWriter(indexFile)));
		} catch (IOException e) {
			U.p("ERROR: Could not create file writer for our report");
			e.printStackTrace();
		}
		
		/* definitely using the whole sequence */
		Properties.useSequenceRegion = false;
		Properties.useOnlyForwardsFrames = false;
		
		/* keep all matches, even bad ones */
		Properties.minimumScore = 1;
		
		/* these properties are fixed for our text sets */
		Properties.numberOfMissedCleavages = 2;
		Properties.precursorTolerance = 2000;
		Properties.fragmentTolerance = 300;
		
		/* this needs to happen or the text reports at the end crash */
		Properties.isSequenceFileDNA = !Properties.testSequenceIsProtein;
		
		/* What scoring mechanism? */
		Properties.scoringMethodName = "Peppy.Match_IMP";
		Properties.matchConstructor = new MatchConstructor(Properties.scoringMethodName);
			
		
		/* save our properties */
		Properties.generatePropertiesFile(reportFolder);
		
	}
	
	
	private static void addTests() {
		//set up which tests we will perform
		String testDirectoryName = Properties.testDirectory.getAbsolutePath();
		tests = new ArrayList<TestSet>();
//		tests.add(new TestSet(testDirectoryName, "ecoli"));
		tests.add(new TestSet(testDirectoryName, "human"));
		tests.add(new TestSet(testDirectoryName, "aurum"));	
		tests.add(new TestSet(testDirectoryName, "USP top 10"));
	}
	
	
	/**
	 * Completes a search on all test sets in our list using the correct (forward) digestion of our protein database.
	 */
	public static void searchTestSets() {
		Sequence sequence;
		if (Properties.testSequenceIsProtein) {
			sequence = new Sequence_Protein(Properties.testSequence);	
		} else {
			sequence = new Sequence_DNA(Properties.testSequence);
		}
		
		/* where we store our peptide chunk */
		ArrayList<Peptide> peptides = sequence.extractMorePeptides(false);
		
		/* find total peptides */
		int totalPeptides = 0;
			
		/* repeat extracting peptides until end of sequence has been reached */
		while (peptides != null) {
			/* sort the peptides */
			Collections.sort(peptides);
			
			/* track database size */
			forwardsDatabaseSize += peptides.size();
			
			/* report on peptides size */
			U.p("peptide count for this batch is: " + peptides.size());
			totalPeptides += peptides.size();
			
			/* get matches for each of our sets for each of these peptides */
			for (TestSet test: tests) {
				test.findPositiveMatches(peptides);
			}
			
			/* get our peptides and report size */
			peptides = sequence.extractMorePeptides(false);
			
		}
		
		U.p("found this many total peptides: " + totalPeptides);
		
		for (TestSet test: tests) {
			test.calculateStastics();
		}
		
	}
	
	
	private static void createListOfPeptidesNotFoundInTheDatabase() {
		U.p("creating list of unfound peptides...");
		
//		Sequence_Protein sequence = new Sequence_Protein(databaseFile);	
		Sequence_DNA sequence = new Sequence_DNA(new File("/Users/risk2/PeppyData/sequences/ecoli/ecoli.fasta"));	
		
		/* where we store our peptide chunk */
		ArrayList<Peptide> peptides = sequence.extractMorePeptides(false);
		
		/* load the correct peptides so we can find which are not in the database */
		ArrayList<Peptide> correctPeptides = new ArrayList<Peptide>();
		for (TestSet test: tests) {
			correctPeptides.addAll(test.loadCorrectPeptides());
		}
		
		/* repeat extracting peptides until end of sequence has been reached */
		while (peptides.size() > 0) {
			
			/* remove the peptides that we do find so that our list of correct peptides contains only those not found */
			for (int i = 0; i < correctPeptides.size(); i++) {
				Peptide peptide = correctPeptides.get(i);
				if (isPeptidePresentInList(peptide, peptides) != -1) {
					correctPeptides.remove(i);
					i--;
					break;
				}
			}
			
			/* get our peptides and report size */
			peptides.clear();
			System.gc();
			peptides = sequence.extractMorePeptides(false);
		}
		
		
		
		File unfoundPeptides = new File(reportFolder,"unfound peptides.txt");
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(unfoundPeptides)));
			/* print the list of unfound peptides */
			for (Peptide peptide: correctPeptides) {
				pw.println(peptide.getAcidSequenceString());
			}
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		

		
	}
	
	
	
	public static void createReport() {
		U.p("Data collected, now generating the report...");

		try {
			/* save raw matches */
			for (TestSet testSet: tests) {
				File testDirectory = new File(reportFolder, testSet.getName());
				testDirectory.mkdirs();
				File matchesFile = new File(testDirectory, "/matches.txt");
				PrintWriter matchPW = new PrintWriter(new BufferedWriter(new FileWriter(matchesFile)));
				ArrayList<Match> matches = testSet.getPositiveMatches();
				Collections.sort(matches);
				for (Match match: matches) {
					matchPW.println(match);
				}
				matchPW.flush();
				matchPW.close();
			}
			
			File indexFile = new File(reportFolder, "index.html");
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(indexFile)));
			NumberFormat nfPercent = NumberFormat.getPercentInstance();
			nfPercent.setMaximumFractionDigits(2);
			
			pw.println("<html>");
			pw.println("<body>");
			pw.println("<h1>Validation Report</h1>");
			
			pw.println("<br>");
			String scoringMethod = Properties.scoringMethodName;
			pw.println("Scoring method: " + scoringMethod);
			
			pw.println("<br>");
			pw.println("Database: " + databaseFile.getName());
			
			pw.println("<h2>Basic performance metrics</h2>");
			pw.println("<table border=1>");
			
			
			//headers
			pw.println("<tr>");
			pw.println("<td>Metric</td>");
			for (TestSet testSet: tests) {
				pw.println("<th>" + testSet.getName() + "</th>");
			}
			
			//# spectra in the set
			pw.println("<tr>");
			pw.println("<td># spectra in the set</td>");
			for (TestSet testSet: tests) {
				pw.println("<td>" + testSet.getSetSize() + "</td>");
			}
			
			//average number of peaks per spectrum
			pw.println("<tr>");
			pw.println("<td>Average # of peaks / spectrum</td>");
			for (TestSet testSet: tests) {
				pw.println("<td>" + testSet.getAverageNumberOfPeaksPerSpectrum() + "</td>");
			}
			
			//Milliseconds per spectrum
			pw.println("<tr>");
			pw.println("<td>Milliseconds per spectrum</td>");
			for (TestSet testSet: tests) {
				pw.println("<td>" + testSet.getMilisecondsPerSpectrum() + "</td>");
			}
			
			
			//# of correct TPs
			pw.println("<tr>");
			pw.println("<td># of true positives</td>");
			for (TestSet testSet: tests) {
				pw.println("<td>" + testSet.getTrueTally() + "</td>");
			}
			
			//% of correct TPs
			pw.println("<tr>");
			pw.println("<td>% of true positives</td>");
			for (TestSet testSet: tests) {
				pw.println("<td>" + (nfPercent.format((double) testSet.getTrueTally() / testSet.getSetSize())) + "</td>");
			}
			
			//# found at 1% FPR
			pw.println("<tr>");
			pw.println("<td># found at 5% FPR</td>");
			for (TestSet testSet: tests) {
				pw.println("<td>" + testSet.getTrueTallyAtFivePercentError() + 
						" (" + nfPercent.format(testSet.getPercentAtFivePercentError()) + ")</td>");
			}
			
//			//% of total found at 1% FPR
//			pw.println("<tr>");
//			pw.println("<td>% of total found at 1% FPR</td>");
//			for (TestSet testSet: tests) {
//				pw.println("<td>" + nfPercent.format(testSet.getPercentAtFivePercentError()) + "</td>");
//			}
//			//E value at 1% FPR
//			pw.println("<tr>");
//			pw.println("<td>E value at 1% FPR</td>");
//			for (TestSet testSet: tests) {
//				pw.println("<td>" + testSet.getEValueAtFivePercentError() + "</td>");
//			}
			
			//PR Curve
			pw.println("<tr>");
			pw.println("<td>Precision-recall curve</td>");
			for (TestSet testSet: tests) {
				BufferedImage bufferedImage = testSet.generatePrecisionRecallCurve();
				File testDirectory = new File(reportFolder, testSet.getName());
				testDirectory.mkdirs();
				File imageFile = new File(testDirectory, "/precision-recall.jpg");
				String fileNameForPRCurve = testSet.getName() + "/" + imageFile.getName();
				ImageIO.write(bufferedImage,"JPG",imageFile);
				pw.println("<td><img src=\"" + fileNameForPRCurve + "\" width=200></td>");
			}
			
			//Area under PR Curve
			pw.println("<tr>");
			pw.println("<td>Area under PR Curve</td>");
			for (TestSet testSet: tests) {
				pw.println("<td>" + nfPercent.format(testSet.getAreaUnderPRCurve()) + "</td>");
			}
			
			pw.println("</table>");		

			

			pw.println("</body>");
			pw.println("</html>");
			
			pw.flush();
			pw.close();
			
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public static void createResultsFiles() {
		for (TestSet test: tests) {
			try {
				File folder = new File(reportFolder, test.getName());
				folder.mkdir();
				File indexFile = new File(folder, test.getName() + " results.txt");
				PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(indexFile)));
				ArrayList<MatchContainer> topForwardsTestedMatches = test.getTestedMatches();
				U.p("number of topForwardsTestedMatches: " + topForwardsTestedMatches.size());
				Match.setSortParameter(Match.SORT_BY_SCORE);
				Collections.sort(topForwardsTestedMatches);
				pw.println("id\tveracity\tscore\tsequence\ttrueScore\ttrueSequence");
				for (MatchContainer match: topForwardsTestedMatches) {
					StringBuffer sb = new StringBuffer();
					sb.append(match.getMatch().getSpectrum().getId());
					sb.append('\t');
					sb.append(match.isTrue());
					sb.append('\t');
					sb.append(match.getMatch().getScore());
					sb.append('\t');
					sb.append(match.getMatch().getPeptide().getAcidSequenceString());
					sb.append('\t');
					sb.append(match.getTrueMatch().getScore());
					sb.append('\t');
					sb.append(match.getTrueMatch().getPeptide().getAcidSequenceString());
					pw.println(sb.toString());
				}
				pw.flush();
				pw.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	
	


	public static void generateRightReport(TestSet test) {
		String testName = test.getName();
		ArrayList<MatchContainer> testedMatches = test.getTestedMatches();
		File testFileFolder = new File(reportFolder, testName);
		testFileFolder.mkdirs();
		File testFile = new File(testFileFolder, "right.html");
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(testFile)));
			
			//print the header of the file
			pw.println("<html><body><table>");
			pw.println("<tr><th>acid</th><th>e value</th><th>score</th>");
			
			for (MatchContainer matchContainer: testedMatches) {
				if (matchContainer.isTrue()) {
					Match ourMatch = matchContainer.getMatch();
					
					pw.println("<tr>");
					
					pw.println("<td>");
					pw.println(ourMatch.getPeptide().getAcidSequenceString());
					pw.println("</td>");
					
					pw.println("<td>");
					pw.println(ourMatch.getScore());
					pw.println("</td>");
					
					pw.println("</tr>");
					pw.println();
					
				}
			}
			
			//print the footer and close out
			pw.println("</table></body></html>");
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	

	
	public static int isPeptidePresentInList(Peptide peptide, ArrayList<Peptide> peptides) {
		/* find the start location */
		int startIndex = MathFunctions.findFirstIndexGreater(peptides, peptide.getMass() - .01);
		startIndex -= 8;
		if (startIndex < 0) startIndex = 0;	
		
		for (int i = startIndex; i < peptides.size(); i++) {
			if (peptide.equals(peptides.get(i))) {
				return i;
			}
			
			/* break out if we are now in the heavier range */
			if (peptides.get(i).getMass() > peptide.getMass()) {
				break;
			}
		}
		return -1;
	}

}
