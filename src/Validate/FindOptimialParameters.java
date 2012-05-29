package Validate;

import java.awt.geom.Point2D;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Random;

import Graphs.PRCurve;
import Peppy.Match;
import Peppy.MatchConstructor;
import Peppy.Peppy;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.Sequence;
import Peppy.Sequence_DNA;
import Peppy.Sequence_Protein;
import Peppy.Spectrum;
import Peppy.MatchesSpectrum;
import Peppy.U;

/**
 * IMPORTANT NOTE:  When running this, be sure to turn all peak cleaning of spectra off.
 * 
 * This walks through precursor tolerance and fragment tolerance to find optimal settings.
 * One thing to keep in mind is that database size affects the optimal precursor tolerance.
 * Usually, the greater the database size, the smaller the optimal tolerance.
 * 
 * @author Brian Risk
 *
 */
public class FindOptimialParameters {
	
	public static void main(String[] args) {	
		findOptimalParameters();
		U.p("done.");
	}
	
	/**
	 * When we don't know what the proper value for the fragment tolerance or precursor tolerance
	 */
	public static void findOptimalParameters() {
		
		//how many missed cleavages when we digest
		Properties.numberOfMissedCleavages = 1;
		
		//What scoring mechanism?
		Properties.scoringMethodName = "Peppy.Match_IMP";
		Properties.matchConstructor = new MatchConstructor(Properties.scoringMethodName);
		
		//set up our tests
		String testDirectoryName = "/Users/risk2/PeppyData/tests/";
//		TestSet test = new TestSet(testDirectoryName, "USP top 10", Color.DARK_GRAY);
		TestSet test = new TestSet(testDirectoryName, "aurum");
		
		//get our peptides
		File databaseFile = new File("/Users/risk2/PeppyData/tests/databases/uniprot_sprot.fasta");
		Sequence_Protein sequence = new Sequence_Protein(databaseFile);	
		ArrayList<Peptide> peptides = sequence.extractAllPeptides(false);
		
		//report thing
		NumberFormat numberFormat = NumberFormat.getInstance();
		numberFormat.setMaximumFractionDigits(2);
		
		//where we are saving it
		File parentDirectory = new File("optimal parameters/" + test.getName());
		parentDirectory.mkdirs();
		
		try {
			PrintWriter pw = new PrintWriter(new FileWriter(new File (parentDirectory, "optimal-parameters-report-005.txt")));
			PrintWriter prGrid = new PrintWriter(new FileWriter(new File (parentDirectory, "PR-grid-00.html5")));
			PrintWriter fprGrid = new PrintWriter(new FileWriter(new File (parentDirectory, "FPR-grid-005.html")));
			prGrid.println("<html><body><table>");
			fprGrid.println("<html><body><table>");
			for (double precursorTolerance = 5; precursorTolerance < 600; precursorTolerance += 5){
				prGrid.println("<tr>");
				fprGrid.println("<tr>");
				for (double fragmentTolerance = 5; fragmentTolerance < 600; fragmentTolerance += 5){
					Properties.precursorTolerance = precursorTolerance;
//					double fragmentTolerance = 0.34;
					Properties.fragmentTolerance = fragmentTolerance;
					test.resetTest();
					test.findPositiveMatches(peptides);
					test.calculateStastics();
					String reportString = 
						numberFormat.format(precursorTolerance) + "," +
						numberFormat.format(fragmentTolerance) + "," +
						test.getTrueTally() + "," +
						test.getPercentAtFivePercentError() + "," +
						test.getAreaUnderPRCurve();
					U.p(reportString);
					pw.println(reportString);
					
					//print the cell in our grids
					prGrid.println("\t<td bgcolor=\"#" + U.getRGBStringFromPercent(test.getAreaUnderPRCurve()) + "\">&nbsp</td>");
					fprGrid.println("\t<td bgcolor=\"#" + U.getRGBStringFromPercent(test.getPercentAtFivePercentError()) + "\">&nbsp</td>");
				}
				prGrid.println("</tr>");
				fprGrid.println("</tr>");
			}
			prGrid.println("</table></body></html>");
			fprGrid.println("</table></body></html>");
			pw.flush();
			pw.close();
			prGrid.flush();
			prGrid.close();
			fprGrid.flush();
			fprGrid.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		

		peptides.clear();
		System.gc();
		
	}
	
	
//	public static void  fop2(String[] args) {
//
//		Peppy.init(args);
//
//		Properties.fragmentTolerance = 100;
//		Properties.precursorTolerance = 100;
//		Properties.minimumScore = 15;
//
//		/* setting objects to format percentages */
//		NumberFormat nfPercent = NumberFormat.getPercentInstance();
//		nfPercent.setMaximumFractionDigits(2);
//		
//		
//		/* Loading a subset of our spectra */
//		ArrayList<File> spectraFiles = new ArrayList<File>();
//		Spectrum.loadSpectraFilesFromFolder(Properties.spectraDirectoryOrFile, spectraFiles);
//		int setSize = Properties.numberOfSpectraToUseForFDR;
//		if (setSize > spectraFiles.size()) setSize = spectraFiles.size();
//		ArrayList<Spectrum> spectra = new ArrayList<Spectrum>();
//		Random random = new Random();
//		File spectrumFile;
//		while (spectra.size() < setSize) {
//			spectrumFile = spectraFiles.remove(random.nextInt(spectraFiles.size()));
//			spectra.addAll(Spectrum.loadSpectra(spectrumFile));
//		}
//		for (int i = 0; i < spectra.size(); i++) {
//			spectra.get(i).setId(i);
//		}
//	
//		
//
//		/* construct a list of points on the PR curve */
//		ArrayList<Point2D.Double> points = new ArrayList<Point2D.Double>();
//		int truePositiveCount = 0;
//		int falsePositiveCount = 0;
//		int decoyAssignmentCount = 0;
//		int targetAssignmentCount = 0;
//		double precision, recall;
//		for (int i = 0; i < allMatches.size(); i++) {
//			Match match = allMatches.get(i);
//			if (match.isDecoy()) {
//				decoyAssignmentCount++;
//			} else {
//				targetAssignmentCount++;
//			}
//			falsePositiveCount = 2 * decoyAssignmentCount;
//			truePositiveCount = (i + 1) - falsePositiveCount;
//			
//			precision = (double) truePositiveCount / (truePositiveCount + falsePositiveCount);
//			recall = (double) truePositiveCount / setSize;
//			Point2D.Double point = new Point2D.Double(recall, precision);
//			points.add(point);
//		}
//		
//
//		
//		/* set a folder to store our reports */
//		String identifierString = Properties.precursorTolerance + "-" + Properties.fragmentTolerance;
//		File reportDir = new File("FDR/" + identifierString + " " + U.getDateString());
//		reportDir.mkdirs();
//		
//		/* save our properties */
//		Properties.generatePropertiesFile(reportDir);
//		
//		/* save curve image */
//		PRCurve prCurve = new PRCurve(points);
//		prCurve.writeFile(new File(reportDir, "pr-" + identifierString + ".jpg"));
//		
//		/* save FDR report */
//		double areaUndrePRCurve = -1;
//		try {
//			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(reportDir, "FDR-" + identifierString + ".txt"))));
//			/* parameters */
//			pw.println("sequenceDirectoryOrFile: " +Properties.sequenceDirectoryOrFile.getAbsolutePath());
//			pw.println("spectraDirectoryOrFile: " +Properties.spectraDirectoryOrFile.getAbsolutePath());
//			pw.println("precursorTolerance: " +Properties.precursorTolerance);
//			pw.println("fragmentTolerance: " +Properties.fragmentTolerance);
//			pw.println();
//			
//			/*print nice stats */
//			areaUndrePRCurve = prCurve.calculateAreaUnderCurve();
//			pw.println("Area under PR curve: " + areaUndrePRCurve);
//			pw.println();
//			
//			/*print first line*/
//			pw.println("Precision\tPercent Found\tIMP value");
//			
//			/* print the rest of the lines */
//			double percent = 1;
//			double increment = 0.01;
//			for (int i = 1; i < points.size(); i++) {
//				if (points.get(i).y < percent && points.get(i - 1).x != points.get(i).x) {
//					percent -= increment;
//					String percentFoundString = (points.get(i - 1).y * 100) + "%";
//					pw.println(percentFoundString + "\t" + points.get(i - 1).x + "\t" + allMatches.get(i - 1).getScore());
//				}
//			}
//			/* to print out the final one */
//			String percentFoundString = (points.get(points.size()  - 1).y * 100) + "%";
//			pw.println(percentFoundString + "\t" + points.get(points.size() - 1).x + "\t" + allMatches.get(points.size() - 1).getScore());
//			
//			pw.flush();
//			pw.close();
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//		
//
//	}
	
	private static ArrayList<Peptide> getTopPeptides(ArrayList<Spectrum> spectra) {
		//Get references to our sequence files -- no nucleotide data is loaded at this point
		ArrayList<Sequence> sequences = Sequence.loadSequenceFiles(Properties.sequenceDirectoryOrFile);
		
		
		/* set up where we will hold all of the matches for our spectra */
		ArrayList<MatchesSpectrum> spectraMatches = new ArrayList<MatchesSpectrum>(spectra.size());
		for (Spectrum spectrum: spectra) {
			spectraMatches.add(new MatchesSpectrum(spectrum));
		}
				
		//getting forwards matches
		U.p("getting forwards matches");
		Peppy.getMatches(sequences, spectraMatches);

		
		//need to initialize things now that we have found matches
		sequences = Sequence_DNA.loadSequenceFiles(Properties.sequenceDirectoryOrFile);
		
		//getting reverse matches -- need to reload the sequences
		U.p("getting reverse matches");
		ArrayList<Match> reverseMatches = Peppy.getDecoyMatches(sequences, spectraMatches);
		reverseMatches = Peppy.reduceMatchesToOnePerSpectrum(reverseMatches);
		Collections.sort(reverseMatches);
		
		/* reducing our peptides to eliminate repeats */
		Hashtable <String, Peptide> hash = new Hashtable<String, Peptide>();
		for (MatchesSpectrum sm: spectraMatches) {
			ArrayList<Match> matches = sm.getMatches();
			for (Match match: matches) {
				Peptide peptide = match.getPeptide();
				hash.put(peptide.getAcidSequenceString(), peptide);
			}
		}
		
		ArrayList<Peptide> out = new ArrayList<Peptide>(hash.values());
		return out;
	}

}
