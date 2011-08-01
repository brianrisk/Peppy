package Validate;

import java.awt.Color;
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
import Utilities.U;

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
		
		Properties.maximumNumberOfMatchesForASpectrum = 1;
		
		//What scoring mechanism?
		Properties.scoringMethodName = "Peppy.Match_IMP";
		Properties.matchConstructor = new MatchConstructor(Properties.scoringMethodName);
		
		//set up our tests
		String testDirectoryName = "/Users/risk2/PeppyOverflow/tests/";
		TestSet test = new TestSet(testDirectoryName, "USP top 10", Color.DARK_GRAY);
		
		//get our peptides
		File databaseFile = new File("/Users/risk2/PeppyOverflow/tests/databases/uniprot_sprot.fasta");
		Sequence_Protein sequence = new Sequence_Protein(databaseFile);	
		ArrayList<Peptide> peptides = sequence.extractAllPeptides(false);
		
		//report thing
		NumberFormat numberFormat = NumberFormat.getInstance();
		numberFormat.setMaximumFractionDigits(2);
		
		try {
			PrintWriter pw = new PrintWriter(new FileWriter(new File ("optimal-parameters-report-005.txt")));
			PrintWriter prGrid = new PrintWriter(new FileWriter(new File ("PR-grid-00.html5")));
			PrintWriter fprGrid = new PrintWriter(new FileWriter(new File ("FPR-grid-005.html")));
			prGrid.println("<html><body><table>");
			fprGrid.println("<html><body><table>");
			for (double precursorTolerance = 0.005; precursorTolerance < 0.6; precursorTolerance += 0.005){
				prGrid.println("<tr>");
				fprGrid.println("<tr>");
				for (double fragmentTolerance = 0.06; fragmentTolerance < .61; fragmentTolerance += 0.01){
					Properties.precursorTolerance = precursorTolerance;
//					double fragmentTolerance = 0.34;
					Properties.fragmentTolerance = fragmentTolerance;
					test.resetTest();
					test.findPositiveMatches(peptides);
					test.cleanMatches();
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

}
