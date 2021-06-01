package Validate;

import Peppy.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

/**
 * IMPORTANT NOTE:  When running this, be sure to turn all peak cleaning of spectra off.
 * <p>
 * This walks through precursor tolerance and fragment tolerance to find optimal settings.
 * One thing to keep in mind is that database size affects the optimal precursor tolerance.
 * Usually, the greater the database size, the smaller the optimal tolerance.
 *
 * @author Brian Risk
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
        TestSet test = new TestSet(testDirectoryName, "aurum");

        //get our peptides
        File databaseFile = new File("/Users/risk2/PeppyData/tests/databases/uniprot_sprot.fasta");
        SequenceProtein sequence = new SequenceProtein(databaseFile);
        ArrayList<Peptide> peptides = sequence.extractAllPeptides(false);

        //report thing
        NumberFormat numberFormat = NumberFormat.getInstance();
        numberFormat.setMaximumFractionDigits(2);

        //where we are saving it
        File parentDirectory = new File("optimal parameters/" + test.getName());
        parentDirectory.mkdirs();

        try {
            PrintWriter pw = new PrintWriter(new FileWriter(new File(parentDirectory, "optimal-parameters-report-005.txt")));
            PrintWriter prGrid = new PrintWriter(new FileWriter(new File(parentDirectory, "PR-grid-00.html5")));
            PrintWriter fprGrid = new PrintWriter(new FileWriter(new File(parentDirectory, "FPR-grid-005.html")));
            prGrid.println("<html><body><table>");
            fprGrid.println("<html><body><table>");
            for (double precursorTolerance = 5; precursorTolerance < 600; precursorTolerance += 5) {
                prGrid.println("<tr>");
                fprGrid.println("<tr>");
                for (double fragmentTolerance = 5; fragmentTolerance < 600; fragmentTolerance += 5) {
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


    private static ArrayList<Peptide> getTopPeptides(ArrayList<Spectrum> spectra) {
        //Get references to our sequence files -- no nucleotide data is loaded at this point
        ArrayList<Sequence> sequences = Sequence.loadSequenceFiles(Properties.sequenceDirectoryOrFile);

        /* set up where we will hold all of the matches for our spectra */
        ArrayList<MatchesSpectrum> spectraMatches = new ArrayList<MatchesSpectrum>(spectra.size());
        for (Spectrum spectrum : spectra) {
            spectraMatches.add(new MatchesSpectrum(spectrum));
        }

        //getting forwards matches
        U.p("getting forwards matches");
        Peppy.getMatches(sequences, spectraMatches);

        //need to initialize things now that we have found matches
        sequences = SequenceNucleotide.loadSequenceFiles(Properties.sequenceDirectoryOrFile);

        //getting reverse matches -- need to reload the sequences
        U.p("getting reverse matches");
        ArrayList<Match> reverseMatches = Peppy.getDecoyMatches(sequences, spectraMatches);
        reverseMatches = Peppy.reduceMatchesToOnePerSpectrum(reverseMatches);
        Collections.sort(reverseMatches);

        /* reducing our peptides to eliminate repeats */
        Hashtable<String, Peptide> hash = new Hashtable<String, Peptide>();
        for (MatchesSpectrum sm : spectraMatches) {
            ArrayList<Match> matches = sm.getMatches();
            for (Match match : matches) {
                Peptide peptide = match.getPeptide();
                hash.put(peptide.getAcidSequenceString(), peptide);
            }
        }

        return new ArrayList<Peptide>(hash.values());
    }

}
