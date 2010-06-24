package Validate;

import java.awt.Color;
import java.awt.Graphics2D;
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

import Peppy.JavaGFS;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.ProteinDigestion;
import Peppy.ScoringThread;
import Peppy.Sequence;
import Peppy.Spectrum;
import Peppy.Match;
import Utilities.U;

public class GenerateValidationReport {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		//Hello, world!
		U.p("Are you ready for the food ball?  I mean: football.  I mean:  validation report");
		
		File validationFolder = new File("Validation");
		validationFolder.mkdir();
		
		File indexFile = new File(validationFolder, "index.html");
		
		//We only want one match per spectrum
		Properties.maximumNumberOfMatchesForASpectrum = 5;
		
		//how many missed cleavages when we digest
		Properties.numberOfMissedCleavages = 2;
		
		//load the swis prot protein database
//		ArrayList<Peptide> peptides = ProteinDigestion.getPeptidesFromProteinFile(new File("/Users/risk2/PeppyOverflow/tests/databases/uniprot_sprot20100614.fasta"));
		ArrayList<Peptide> peptides = ProteinDigestion.getPeptidesFromProteinFile(new File("/Users/risk2/PeppyOverflow/tests/databases/uniprot_sprot.fasta"));
//		Sequence ecoliSequence = new Sequence(new File("/Users/risk2/PeppyOverflow/sequences ecoli/ecoli.fasta"));
//		ArrayList<Peptide> peptides = ecoliSequence.extractPeptides();
		
		//set up which tests we will perform
		ArrayList<String> tests = new ArrayList<String>();
//		tests.add("human");
//		tests.add("ecoli");
//		tests.add("aurum");
		tests.add("USP");
		
		
		//get the matches for each of our tests
		for (int i = 0; i < tests.size(); i++) {
			String test = tests.get(i);
			U.p("Getting matches for: " + test);
			
			//load spectra for this test
			ArrayList<Spectrum> spectra = Spectrum.loadSpectraFromFolder("/Users/risk2/PeppyOverflow/tests/" + test + "/spectra");
			U.p("loaded " +spectra.size() + " spectra.");
			
			//get the matches
			ArrayList<Match> matches = JavaGFS.asynchronousDigestion(peptides, spectra, null);
			
			//See which of these matches are true
			ArrayList<MatchContainer> testedMatches = new ArrayList<MatchContainer>(matches.size());
			for (Match match: matches) {
				if (match == null) U.p("null here");
				testedMatches.add(new MatchContainer(match));
			}
			
			
			//Sort by e value
			Collections.sort(testedMatches);
			int trueTally = 0;
			int falseTally = 0;
			boolean onePercenThresholdHasBeenReached = false;
			boolean firstFalseFound = false;
			double firstFalseEValue = 0;
			for (MatchContainer match: testedMatches) {
				//U.p(match.isTrue() + " " + match.getEValue());
				if (match.isTrue()) {
					trueTally++;
				} else {
					falseTally++;
					if (!firstFalseFound) {
						firstFalseFound = true;
						firstFalseEValue = match.getEValue();
						U.p("First fasle found at this e Value: " + firstFalseEValue);
					}
				}
				if (!onePercenThresholdHasBeenReached) {
					if ((double) falseTally / testedMatches.size() >= 0.01) {
						onePercenThresholdHasBeenReached = true;
						U.p("Number at one percent error rate: " + trueTally);
						U.p("E value at this rate is: " + match.getEValue());
					}
				}
			}
			U.p("True hit tally: " + trueTally);
			U.p("True hit percent: " + (double) trueTally / testedMatches.size());
			U.p();
			
			generateWrongReport(test, testedMatches, peptides);
			
			
//			generateReport(test, testedMatches);
		}
		 
		//get all the matches		

	}
	
	public static void generateWrongReport(String test, ArrayList<MatchContainer> testedMatches, ArrayList<Peptide> peptides) {
		//intro!
		U.p("printing out all the spectra that we got wrong and comparing the two different matches");
		
		int width = 1000;
		int height = 200;
		
		File testFileFolder = new File(Properties.reportDirectory, test);
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
					try {
					correctPeptide = new Peptide(matchContainer.getCorrectAcidSequence());
					} catch (Exception e) {
						U.p("Exception!");
						U.p("is acidSequence null? answer: " + matchContainer.getCorrectAcidSequence() == null);
						U.p("spectrum file: " + matchContainer.getMatch().getSpectrum().getFile().getName());
					}
					Match ourMatch = matchContainer.getMatch();
					Match trueMatch = matchContainer.getTrueMatch();
					pw.println("<p>");
					pw.println("E value: " + ourMatch.getEValue());
					pw.println("<br>");
					pw.println("Spectrum file name: " + ourMatch.getSpectrum().getFile().getName());
					pw.println("<br>");
					pw.println("Correct peptide is in the database: " + isPeptidePresentInList(correctPeptide,peptides));
					pw.println("<br>");
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
					
					//print int out
					pw.println("ours vs theirs: " + ourMatch.getPeptide().getAcidSequence() + " vs " + matchContainer.getCorrectAcidSequence() + "<br>");
					pw.print("<img src=\"");
					pw.print(imagesFolder.getName() + "/" + ourMatchFile.getName());
					pw.print("\"><br>");
					pw.print("<img src=\"");
					pw.print(imagesFolder.getName() + "/" + trueMatchFile.getName());
					pw.print("\"><br>");
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
	
	public static void generateReport(String test, ArrayList<MatchContainer> testedMatches) {
		File folder = new File ("validaiont/" + test);
		folder.mkdirs();
		
	}
	
	public static void generatePrecisionRecallCurve(ArrayList<MatchContainer> testedMatches, int width, int height) {
		int size = testedMatches.size();
		int trueCount = 0;
		double previousPrecision = 0;
		double previousRecall = 0;
		double precision = 0; 
		double recall = 0;
		
		int xLoc, yLoc;
		
		//setting up Graphics context
		BufferedImage bdest = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		Graphics2D g = bdest.createGraphics();
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_SPEED));
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF));
		g.setColor(Color.white);
		g.fillRect(0,0,width,height);
		
		
		for(int i = 0; i < size; i++) {
			MatchContainer match = testedMatches.get(i);
			if (match.isTrue()) {
				trueCount++;
			}
			previousPrecision = precision;
			previousRecall = recall;
			precision = (double) trueCount / (i + 1);
			recall = (double) trueCount / size;	
			
		}
		
		try {
			ImageIO.write(bdest,"JPG",new File("out.jpg"));
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static boolean isPeptidePresentInList(Peptide peptide, ArrayList<Peptide> peptides) {
		int peptideIndex = ScoringThread.findFirstIndexWithGreaterMass(peptides, peptide.getMass() - .01);
		double peptideMassButBigger = peptide.getMass() + .01;
		for (int i = peptideIndex; i < peptides.size(); i++) {
			if (peptide.getAcidSequence().equals(peptides.get(i).getAcidSequence())) {
				return true;
			}
			if (peptides.get(i).getMass() > peptideMassButBigger) {
				break;
			}
		}
		return false;
	}

}
