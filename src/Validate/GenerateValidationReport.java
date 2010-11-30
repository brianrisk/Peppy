package Validate;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import Peppy.Match;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.ProteinDigestion;
import Peppy.ScoringThread;
import Utilities.U;

public class GenerateValidationReport {
	
	public static ArrayList<TestSet> tests;
	public static File databaseFile;
	public static PrintWriter indexWriter;
	public static boolean doReverse = false;

	/**
	 * @param args
	 */
	public static void main(String[] args) {	
		setUp();
		forwards();
		if (doReverse) reverse();
		createReport();
		U.p("done.");
	}
	
	public static void setUp() {
		//Hello, world!
		U.p("Are you ready for the food ball?  I mean: football.  I mean:  validation report");
		Properties.validationDirectory.mkdirs();
		File indexFile = new File(Properties.validationDirectory, "index.html");
		try {
			indexWriter = new PrintWriter(new BufferedWriter(new FileWriter(indexFile)));
		} catch (IOException e) {
			U.p("ERROR: Could not create file writer for our report");
			e.printStackTrace();
		}
		
//		Properties.peakDifferenceThreshold = 0.3;
		
		//how many missed cleavages when we digest
		Properties.numberOfMissedCleavages = 2;
		
		Properties.spectrumToPeptideMassError = 2.0;
		Properties.peakDifferenceThreshold = 0.5;
		
		//we'd prefer not to have duplicate matches -- especially for the correct ones
		Properties.reduceDuplicateMatches = true;
		
		//What scoring mechanism?
//		Properties.defaultScore = Properties.DEFAULT_SCORE_TANDEM_FIT;
		Properties.defaultScore = Properties.DEFAULT_SCORE_HMM;
		HMMScore.HMMClass.HmmSetUp();
//		Properties.highIntensityCleaning = true;
		Properties.localMaximaCleaning = true;
		
		databaseFile = new File("/Users/risk2/PeppyOverflow/tests/databases/uniprot_sprot.fasta");
//		databaseFile = new File("uniprot_sprot.fasta");
		
		//set up which tests we will perform
		tests = new ArrayList<TestSet>();
		
		String testDirectoryName = "/Users/risk2/PeppyOverflow/tests/";
//		String testDirectoryName = "tests/";
		tests.add(new TestSet(testDirectoryName, "ecoli"));
		tests.add(new TestSet(testDirectoryName, "human"));
		tests.add(new TestSet(testDirectoryName, "aurum"));	
		tests.add(new TestSet(testDirectoryName, "USP"));

	}
	
	
	/**
	 * Completes a search on all test sets in our list using the correct (forward) digestion of our protein database.
	 */
	public static void forwards() {
		Properties.maximumNumberOfMatchesForASpectrum = 10;
		//load the peptides
		ArrayList<Peptide> peptides = ProteinDigestion.getPeptidesFromDatabase(databaseFile);
		U.p("forwards database size: " + peptides.size());
//		Sequence ecoli = new Sequence("/Users/risk2/PeppyOverflow/sequences ecoli/ecoli.fasta");
//		ArrayList<Peptide> peptides = ecoli.extractAllPeptides();

		for (TestSet test: tests) {
			U.p("Getting matches for: " + test.getName());
			test.findPositiveMatches(peptides);
		}	
	}
	
	/**
	 * Get the reverse of the database.  Should produce a database of about the same size but
	 * with, most likely, nearly no correct matches.
	 */
	public static void reverse() {
		//We only want one match per spectrum
		Properties.maximumNumberOfMatchesForASpectrum = 1;
		
		ArrayList<Peptide> peptides = ProteinDigestion.getReversePeptidesFromFASTA(databaseFile);
		
		for (TestSet test: tests) {
			U.p("Getting false matches for: " + test.getName());
			test.findFalsePositiveMatches(peptides);
		}		
	}
	
	public static void createReport() {
		U.p("Data collected, now generating the report...");

		try {
			File indexFile = new File(Properties.validationDirectory, "index.html");
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(indexFile)));
			
			pw.println("<html>");
			pw.println("<body>");
			pw.println("<h1>Validation Report</h1>");
			pw.println("<h2>Basic performance metrics</h2>");
			pw.println("<table border=1>");
			
			
			//headers
			pw.println("<tr>");
			pw.println("<td>Metric</td>");
			for (TestSet testSet: tests) {
				pw.println("<td>" + testSet.getName() + "</td>");
			}
			
			//# spectra in the set
			pw.println("<tr>");
			pw.println("<td># spectra in the set</td>");
			for (TestSet testSet: tests) {
				pw.println("<td>" + testSet.getSetSize() + "</td>");
			}
			
			//Time to complete
			pw.println("<tr>");
			pw.println("<td>Time to complete</td>");
			for (TestSet testSet: tests) {
				pw.println("<td>" + testSet.getTimeToComplete() + "</td>");
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
				pw.println("<td>" + ((double) testSet.getTrueTally() / testSet.getSetSize()) + "</td>");
			}
			
			//# found at 1% FPR
			pw.println("<tr>");
			pw.println("<td># found at 1% FPR</td>");
			for (TestSet testSet: tests) {
				pw.println("<td>" + testSet.getTrueTallyAtOnePercentError() + "</td>");
			}
			
			//% of total found at 1% FPR
			pw.println("<tr>");
			pw.println("<td>% of total found at 1% FPR</td>");
			for (TestSet testSet: tests) {
				pw.println("<td>" + testSet.getPercentAtOnePercentError() + "</td>");
			}
			//E value at 1% FPR
			pw.println("<tr>");
			pw.println("<td>E value at 1% FPR</td>");
			for (TestSet testSet: tests) {
				pw.println("<td>" + testSet.getEValueAtOnePercentError() + "</td>");
			}
			
			//PR Curve
			pw.println("<tr>");
			pw.println("<td>Precision-recall curve</td>");
			for (TestSet testSet: tests) {
				pw.println("<td><img src=\"" + testSet.getFileNameForPRCurve() + "\" width=200></td>");
			}
			
			//Area under PR Curve
			pw.println("<tr>");
			pw.println("<td>Area under PR Curve</td>");
			for (TestSet testSet: tests) {
				pw.println("<td>" + testSet.getAreaUnderPRCurve() + "</td>");
			}
			
			
			
			
			pw.println("</table>");
			pw.println("<h2>E value distribution for forwards and reverse database search</h2>");
			
			pw.println("<h3>Forwards</h3>");
			pw.println("<table border=1>");
			pw.println("<tr>");
			pw.println("<td>E value marking top 5% (forwards)</td>");
			for (TestSet testSet: tests) {
				pw.println("<td>" + testSet.getEValueAtPercentForwards(0.05) + "</td>");
			}
			pw.println("<tr>");
			pw.println("<td>E value marking top 25% (forwards)</td>");
			for (TestSet testSet: tests) {
				pw.println("<td>" + testSet.getEValueAtPercentForwards(0.25) + "</td>");
			}
			pw.println("<tr>");
			pw.println("<td>E value marking top 50% (forwards)</td>");
			for (TestSet testSet: tests) {
				pw.println("<td>" + testSet.getEValueAtPercentForwards(0.50) + "</td>");
			}
			pw.println("<tr>");
			pw.println("<td>E value marking top 75% (forwards)</td>");
			for (TestSet testSet: tests) {
				pw.println("<td>" + testSet.getEValueAtPercentForwards(0.75) + "</td>");
			}
			pw.println("<tr>");
			pw.println("<td>E value marking top 95% (forwards)</td>");
			for (TestSet testSet: tests) {
				pw.println("<td>" + testSet.getEValueAtPercentForwards(0.95) + "</td>");
			}
			pw.println("</table>");
			
			//REVERSE
			if (doReverse) {
				pw.println("<h3>Reverse</h3>");
				pw.println("<table border=1>");
				pw.println("<tr>");
				pw.println("<td>E value marking top 5% (reverse)</td>");
				for (TestSet testSet: tests) {
					pw.println("<td>" + testSet.getEValueAtPercentReverse(0.05) + "</td>");
				}
				pw.println("<tr>");
				pw.println("<td>E value marking top 25% (reverse)</td>");
				for (TestSet testSet: tests) {
					pw.println("<td>" + testSet.getEValueAtPercentReverse(0.25) + "</td>");
				}
				pw.println("<tr>");
				pw.println("<td>E value marking top 50% (reverse)</td>");
				for (TestSet testSet: tests) {
					pw.println("<td>" + testSet.getEValueAtPercentReverse(0.50) + "</td>");
				}
				pw.println("<tr>");
				pw.println("<td>E value marking top 75% (reverse)</td>");
				for (TestSet testSet: tests) {
					pw.println("<td>" + testSet.getEValueAtPercentReverse(0.75) + "</td>");
				}
				pw.println("<tr>");
				pw.println("<td>E value marking top 95% (reverse)</td>");
				for (TestSet testSet: tests) {
					pw.println("<td>" + testSet.getEValueAtPercentReverse(0.95) + "</td>");
				}
				pw.println("</table>");
			}

			pw.println("</table>");
			pw.println("</body>");
			pw.println("</html>");
			
			pw.flush();
			pw.close();
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
//		//Visualize the ion matches for all incorrect matches
//		for (TestSet test: tests) {
//			generateWrongReport(test);
//		}
//		
//		
//		//Report on reverse database
//		for (TestSet test: tests) {
//			U.p();
//			U.p(test.getName());
//			ArrayList<Match> falsePositiveMatches = test.getFalsePositiveMatches();
//			Match.setSortParameter(Match.SORT_BY_E_VALUE);
//			Collections.sort(falsePositiveMatches);
//			int testSize = falsePositiveMatches.size();
//			int level05 = (int) (testSize * 0.05);
//			int level25 = (int) (testSize * 0.25);
//			int level50 = (int) (testSize * 0.50);
//			int level75 = (int) (testSize * 0.75);
//			int level95 = (int) (testSize * 0.95);
//			U.p(" 5%: "+ falsePositiveMatches.get(level05).getEValue());
//			U.p("25%: "+ falsePositiveMatches.get(level25).getEValue());
//			U.p("50%: "+ falsePositiveMatches.get(level50).getEValue());
//			U.p("75%: "+ falsePositiveMatches.get(level75).getEValue());
//			U.p("95%: "+ falsePositiveMatches.get(level95).getEValue());
//		}
	}
	
	
	/**
	 * Used to validate our test set.  We want visualizations of
	 * what the "wrong" peptide aligned with the spectrum compared
	 * to the "right" one.
	 */
	public static void makeOnlyWrongReport() {
		setUp();
		U.startStopwatch();
		Properties.maximumNumberOfMatchesForASpectrum = 1;
		ArrayList<Peptide> peptides = ProteinDigestion.getPeptidesFromFASTA(databaseFile);
		for (TestSet test: tests) {
			test.findPositiveMatches(peptides);
			generateWrongReport(test);
		}
		U.p();
		U.stopStopwatch();
	}

	public static void generateRightReport(TestSet test) {
		String testName = test.getName();
		ArrayList<MatchContainer> testedMatches = test.getTestedMatches();
		File testFileFolder = new File(Properties.validationDirectory, testName);
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
					pw.println(ourMatch.getPeptide().getAcidSequence());
					pw.println("</td>");
					
					pw.println("<td>");
					pw.println(ourMatch.getEValue());
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
	
	public static void generateWrongReport(TestSet test) {
		String testName = test.getName();
		ArrayList<MatchContainer> testedMatches = test.getTestedMatches();
		U.p("printing out all the spectra that we got wrong and comparing the two different matches");
		
		int width = 1000;
		int height = 200;
		
		File testFileFolder = new File(Properties.validationDirectory, testName);
		testFileFolder.mkdirs();
		File testFile = new File(testFileFolder, "wrong.html");
		File imagesFolder = new File(testFileFolder, "wrongImages");
		imagesFolder.mkdir();
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(testFile)));
			
			//print the header of the file
			pw.println("<html><body>");
			
			int imageIndex = 0;
			for (MatchContainer matchContainer: testedMatches) {
				if (!matchContainer.isTrue()) {
					Peptide correctPeptide = null;
					correctPeptide = new Peptide(matchContainer.getCorrectAcidSequence());
					Match ourMatch = matchContainer.getMatch();
					Match trueMatch = matchContainer.getTrueMatch();
					pw.println("<p>");
					pw.println("E value: " + ourMatch.getEValue());
					pw.println("<br>");
					pw.println("Spectrum file name: " + ourMatch.getSpectrum().getFile().getName());
					pw.println("<br>");
//					if (checkToSeeIfPeptideIsInDatabase) {
//						pw.println("Correct peptide is in the database: " + isPeptidePresentInList(correctPeptide, peptides));
//						pw.println("<br>");
//					}
					pw.println("our score: " + ourMatch.ionMatchTally + ", " + ourMatch.getScore());
					pw.println("<br>");
					pw.println("their score: " + trueMatch.ionMatchTally + ", " + trueMatch.getScore());
					pw.println("<br>");
					
					//creating the image of our match
					File ourMatchFile = new File (imagesFolder, imageIndex + "a.jpg");
					SpectralVisualizer.SpectralVisualizer.markMatchingIons(ourMatch.getSpectrum(), ourMatch.getPeptide());
					SpectralVisualizer.SpectralVisualizer.drawSpectrum(ourMatch.getSpectrum(), width, height, ourMatchFile, false);
					
					//creating the image of their "true" match.  Whatever!
					File trueMatchFile = new File (imagesFolder, imageIndex + "b.jpg");
					SpectralVisualizer.SpectralVisualizer.markMatchingIons(ourMatch.getSpectrum(), correctPeptide);
					SpectralVisualizer.SpectralVisualizer.drawSpectrum(ourMatch.getSpectrum(), width, height, trueMatchFile, false);
					
					//don't forget to increment that image index, Brian!
					imageIndex++;
					
					//print in out
					pw.println("ours vs theirs: " + ourMatch.getPeptide().getAcidSequence() + " vs " + matchContainer.getCorrectAcidSequence() + "<br>");
					pw.print("<a href=\"");
					pw.print(imagesFolder.getName() + "/" + ourMatchFile.getName());
					pw.print("\">");
					pw.print("<img src=\"");
					pw.print(imagesFolder.getName() + "/" + ourMatchFile.getName());
					pw.print("\" border=0></a><br>");
					
					pw.print("<a href=\"");
					pw.print(imagesFolder.getName() + "/" + trueMatchFile.getName());
					pw.print("\">");
					pw.print("<img src=\"");
					pw.print(imagesFolder.getName() + "/" + trueMatchFile.getName());
					pw.print("\" border=0></a><br>");
					pw.print("<hr>");
					pw.println();
					pw.println();
				}
			}
			
			//print the footer and close out
			pw.println("</body></html>");
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	

	
	public static int isPeptidePresentInList(Peptide peptide, ArrayList<Peptide> peptides) {
		int peptideIndex = ScoringThread.findFirstIndexWithGreaterMass(peptides, peptide.getMass() - .01);
		double peptideMassButBigger = peptide.getMass() + .01;
		for (int i = peptideIndex; i < peptides.size(); i++) {
			if (peptide.getMass() == peptides.get(i).getMass()) {
				if (peptide.getAcidSequence().equals(peptides.get(i).getAcidSequence())) {
					return i;
				}
			}
			if (peptides.get(i).getMass() > peptideMassButBigger) {
				break;
			}
		}
		return -1;
	}

}
