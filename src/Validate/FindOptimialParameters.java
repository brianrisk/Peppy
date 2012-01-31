package Validate;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.ArrayList;

import Peppy.MatchConstructor;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.Sequence_Protein;
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
//		findOptimalParameters();
		findOptimalTandemFitExponent();
		U.p("done.");
	}
	
	/**
	 * When we don't know what the proper value for the fragment tolerance or precursor tolerance
	 */
	public static void findOptimalParameters() {
		
		//how many missed cleavages when we digest
		Properties.numberOfMissedCleavages = 1;
		
		Properties.maximumNumberOfMatchesForASpectrum = 1;
		
		//What scoring mechanism?
		Properties.scoringMethodName = "Peppy.Match_IMP";
		Properties.matchConstructor = new MatchConstructor(Properties.scoringMethodName);
		
		//set up our tests
		String testDirectoryName = "/Users/risk2/PeppyOverflow/tests/";
//		TestSet test = new TestSet(testDirectoryName, "USP top 10", Color.DARK_GRAY);
		TestSet test = new TestSet(testDirectoryName, "aurum");
		
		//get our peptides
		File databaseFile = new File("/Users/risk2/PeppyOverflow/tests/databases/uniprot_sprot.fasta");
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
	
	
	/**
	 * When we don't know what the proper value for the fragment tolerance or precursor tolerance
	 */
	public static void findOptimalTandemFitExponent() {
		
		//how many missed cleavages when we digest
		Properties.numberOfMissedCleavages = 1;
		
		Properties.maximumNumberOfMatchesForASpectrum = 1;
		
		//What scoring mechanism?
		Properties.scoringMethodName = "Peppy.Match_TandemFit";
		Properties.matchConstructor = new MatchConstructor(Properties.scoringMethodName);
		
		//set up our tests
		String testDirectoryName = "/Users/risk2/PeppyOverflow/tests/";
		TestSet test = new TestSet(testDirectoryName, "human");
		
		//get our peptides
		File databaseFile = new File("/Users/risk2/PeppyOverflow/tests/databases/uniprot_sprot.fasta");
		Sequence_Protein sequence = new Sequence_Protein(databaseFile);	
		ArrayList<Peptide> peptides = sequence.extractAllPeptides(false);
		U.p("peptides created");
		
		//report thing
		NumberFormat numberFormat = NumberFormat.getInstance();
		numberFormat.setMaximumFractionDigits(5);
		
		//where we are saving it
		File parentDirectory = new File("optimal parameters/" + test.getName());
		parentDirectory.mkdirs();
		
		try {
			PrintWriter pw = new PrintWriter(new FileWriter(new File (parentDirectory, "aryaScore-exponent-performance.txt")));
			for (double peakIntensityExponent = 0.00; peakIntensityExponent < 1.0; peakIntensityExponent += 0.01){


				Properties.peakIntensityExponent = peakIntensityExponent;
				test.resetTest();
				test.findPositiveMatches(peptides);
				test.calculateStastics();
				String reportString = 
					numberFormat.format(peakIntensityExponent) + "," +
					test.getTrueTally() + "," +
					test.getPercentAtFivePercentError() + "," +
					test.getAreaUnderPRCurve();
				U.p(reportString);
				pw.println(reportString);

			}

			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		

		peptides.clear();
		System.gc();
		
	}

}
