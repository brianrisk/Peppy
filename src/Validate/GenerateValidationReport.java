package Validate;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Polygon;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;

import javax.imageio.ImageIO;

import Peppy.Match;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.ProteinDigestion;
import Peppy.ScoringThread;
import Utilities.U;

public class GenerateValidationReport {
	
	public static ArrayList<TestSet> tests;
	public static long timeElapsed;
	public static File databaseFile;
	public static PrintWriter indexWriter;


	/**
	 * @param args
	 */
	public static void main(String[] args) {	
		setUp();
		U.startStopwatch();
//		forwards();
//		reverse();
		//Visualize the ion matches for all incorrect matches
		//load the peptides
		Properties.maximumNumberOfMatchesForASpectrum = 1;
		ArrayList<Peptide> peptides = ProteinDigestion.getPeptidesFromProteinFile(databaseFile);
		for (TestSet test: tests) {
			test.findPositiveMatches(peptides);
			generateWrongReport(test);
		}
		U.p();
		timeElapsed = U.stopStopwatch();
//		createReport();
	}
	
	
	public static void setUp() {
		//Hello, world!
		U.p("Are you ready for the food ball?  I mean: football.  I mean:  validation report");
		
		File indexFile = new File(Properties.validationDirectory, "index.html");
		try {
			indexWriter = new PrintWriter(new BufferedWriter(new FileWriter(indexFile)));
		} catch (IOException e) {
			U.p("ERROR: Could not create file writer for our report");
			e.printStackTrace();
		}
		
		//how many missed cleavages when we digest
		Properties.numberOfMissedCleavages = 2;
		
		//we'd prefer not to have duplicate matches -- especially for the correct ones
		Properties.reduceDuplicateMatches = true;
		
		//What scoring mechanism?
		Properties.defaultScore = Properties.DEFAULT_SCORE_TANDEM_FIT;
//		Properties.defaultScore = Properties.DEFAULT_SCORE_HMM;
//		HMMScore.HMMClass.HmmSetUp();
		
		databaseFile = new File("/Users/risk2/PeppyOverflow/tests/databases/uniprot_sprot.fasta");
		
		//set up which tests we will perform
		tests = new ArrayList<TestSet>();

//		U.p("Getting matches for: ecoli");
//		tests.add(new TestSet("ecoli"));
//		
//		U.p("Getting matches for: human");
//		tests.add(new TestSet("human"));
//	
//		U.p("Getting matches for: aurum");
//		tests.add(new TestSet("aurum"));
		
		U.p("Getting matches for: USP");
		tests.add(new TestSet("USP"));
	}
	
	
	/**
	 * Completes a search on all test sets in our list using the correct (forward) digestion of our protein database.
	 */
	public static void forwards() {
		Properties.maximumNumberOfMatchesForASpectrum = 50;
		//load the peptides
		ArrayList<Peptide> peptides = ProteinDigestion.getPeptidesFromProteinFile(databaseFile);

		for (TestSet test: tests) {
			U.p("Getting false matches for: " + test.getName());
			test.findFalsePositiveMatches(peptides);
		}	
	}
	
	
	public static void reverse() {
		//5, 25, 75, 95
		U.p("now for the reverse database...");
		//We only want one match per spectrum
		Properties.maximumNumberOfMatchesForASpectrum = 1;
		
		//Get the reverse of the database.  Should produce a database of about the same size but
		//with, most likely, nearly no correct matches.
		ArrayList<Peptide> peptides = ProteinDigestion.getReversePeptidesFromProteinFile(databaseFile);
		
		for (TestSet test: tests) {
			U.p("Getting false matches for: " + test.getName());
			test.findFalsePositiveMatches(peptides);
		}		
	}
	
	public static void createReport() {
		U.p("Data collected, now generating the report...");
		
		
		//Finding basic stats
		for (TestSet test: tests) {
			int trueTally = 0;
			int falseTally = 0;
			int trueTallyAtOnePercentError = -1;
			double eValueAtOnePercentError = -1;
			ArrayList<MatchContainer> testedMatches = test.getTestedMatches();
			boolean onePercenThresholdHasBeenReached = false;
			for (MatchContainer match: testedMatches) {
				if (match.isTrue()) {
					trueTally++;
				} else {
					falseTally++;
				}
				if (!onePercenThresholdHasBeenReached) {
					if ((double) falseTally / test.getSetSize() >= 0.01) {
						onePercenThresholdHasBeenReached = true;
						trueTallyAtOnePercentError =  trueTally;
						eValueAtOnePercentError = match.getEValue();
					}
				}
			}
			

		}
		
		//Precision recall
		for (TestSet test: tests) {
			generatePrecisionRecallCurve(test, 300, 300);
		}
		
		//Visualize the ion matches for all incorrect matches
		for (TestSet test: tests) {
			generateWrongReport(test);
		}
		
		
		//Report on reverse database
		for (TestSet test: tests) {
			U.p();
			U.p(test.getName());
			ArrayList<Match> falsePositiveMatches = test.getFalsePositiveMatches();
			Match.setSortParameter(Match.SORT_BY_E_VALUE);
			Collections.sort(falsePositiveMatches);
			int testSize = falsePositiveMatches.size();
			int level05 = (int) (testSize * 0.05);
			int level25 = (int) (testSize * 0.25);
			int level50 = (int) (testSize * 0.50);
			int level75 = (int) (testSize * 0.75);
			int level95 = (int) (testSize * 0.95);
			U.p(" 5%: "+ falsePositiveMatches.get(level05).getEValue());
			U.p("25%: "+ falsePositiveMatches.get(level25).getEValue());
			U.p("50%: "+ falsePositiveMatches.get(level50).getEValue());
			U.p("75%: "+ falsePositiveMatches.get(level75).getEValue());
			U.p("95%: "+ falsePositiveMatches.get(level95).getEValue());
		}
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
					pw.println("our score: " + ourMatch.ionMatchTally + ", " + ourMatch.getScoreTandemFit());
					pw.println("<br>");
					pw.println("their score: " + trueMatch.ionMatchTally + ", " + trueMatch.getScoreTandemFit());
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
	
	public static void generatePrecisionRecallCurve(TestSet test, int width, int height) {
		String testName = test.getName();
		ArrayList<MatchContainer>testedMatches = test.getTestedMatches();

		//setting up Graphics context
		BufferedImage bdest = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		Graphics2D g = bdest.createGraphics();
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY));
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
		g.setColor(Color.white);
		g.fillRect(0,0,width,height);
		
		//setting the line color
		g.setColor(Color.red);
		
		int x1 = 0;
		int x2 = 0;
		int y1 = 0;
		int y2 = 0;
		int trueCount = 0;
		double precision = 0; 
		double recall = 0;
		double recallPrevious = 0;
		double area = 0;
		for(int i = 0; i < testedMatches.size(); i++) {
			MatchContainer match = testedMatches.get(i);
			if (match.isTrue()) {
				trueCount++;
			}
			//advancing
			x1 = x2 + 1;
			y1 = y2;
			
			precision = (double) trueCount / (i + 1);
			recallPrevious = recall;
			recall = (double) trueCount / test.getSetSize();	
			
			area += (recall - recallPrevious) * precision;
			
			x2 = (int) (recall * width);
			y2 = (int) ((1.0 - precision) * height);
			
			//in case we are moving so little we are not advancing
			if (x1 <= x2) {
				continue;
			} else {
//				U.p(trueCount + ", " + precision);
				g.setColor(Color.red);
				g.setStroke(new BasicStroke(2.0f));
				g.drawLine(x1, y1, x2, y2);
				//let's fill in the area under the line, yes?
				Polygon polygon = new Polygon();
				polygon.addPoint(x1, y1);
				polygon.addPoint(x2, y2);
				polygon.addPoint(x2, height);
				polygon.addPoint(x1, height);
				g.setColor(new Color(256,0,0,128));
				g.fillPolygon(polygon);
			}
			
		}
		
		U.p("The area is: " + area);
		
//		g.drawLine(x2, y2, width, height);
		
		try {
			ImageIO.write(bdest,"JPG",new File(Properties.validationDirectory, testName + "-precision-recall.jpg"));
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
