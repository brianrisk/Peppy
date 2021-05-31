package Validate;

import Math.MassError;
import Peppy.*;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Test set's have some basic properties which would be cool to know
 * such as:
 * 1) In what ways does the precursor mass differ from the predicted mass
 * 2) are there consistent errors between predicted peaks and observed?
 * <p>
 * If there are consistencies, then the scoring mechanisms can be tuned
 * to accommodate the eccentricities of these results.
 *
 * @author Brian Risk
 */
public class TestSetCharacteristics {

    private String testName;
    private ArrayList<Spectrum> spectra;
    ArrayList<Match> correctMatches;

    public static void main(String args[]) {
        new TestSetCharacteristics("USP");
    }

    public TestSetCharacteristics(String testName) {
        this.testName = testName;

        //load spectra for this test
        spectra = SpectrumLoader.loadSpectra("/Users/risk2/PeppyData/tests/" + testName + "/spectra");

        //set up correct matches
        correctMatches = loadCorrectMatches();

        printIonDifferences(correctMatches);


//		printPrecursorDifferenceToPredictedMass();


        //U.p(countBadPrecursors());

        U.p("done");
    }

    public static void printIonDifferences(ArrayList<Match> matches) {
        double y_ionCount = 0;
        double y_ionDifferenceSum = 0;
        double b_ionCount = 0;
        double b_ionDifferenceSum = 0;
        for (Match match : matches) {
            Peptide peptide = match.getPeptide();
            Spectrum spectrum = match.getSpectrum();
            ArrayList<Peak> peaks = spectrum.getPeaks();
            String peptideString = peptide.getAcidSequenceString();

            //If precursor difference is too much, exit
            double difference = spectrum.getMass() - peptide.getMass();
            if (Math.abs(difference) > 2) continue;

            int i;
            double theoreticalPeakMass, peakMass;


            //we want -1 because most of these spectra will have a match with
            //the last theoretical peak
            int peptideLengthMinusOne = peptideString.length() - 1;

            double[] bIonMatchesWithHighestIntensity = new double[peptideString.length()];
            for (i = 0; i < peptideString.length(); i++) bIonMatchesWithHighestIntensity[i] = 0.0;
            double[] yIonMatchesWithHighestIntensity = new double[peptideString.length()];
            for (i = 0; i < peptideString.length(); i++) yIonMatchesWithHighestIntensity[i] = 0.0;

            double theoreticalPeakLeft;
            double theoreticalPeakRight;
            double peakMatchWithGreatestIntensity = 0;
            double maxIntensity;
            double peakIntensity;

            /* y-ion  */
            theoreticalPeakMass = peptide.getMass() + Properties.rightIonDifference;
            for (i = 0; i < peptideLengthMinusOne; i++) {
                theoreticalPeakMass -= AminoAcids.getWeightMono(peptideString.charAt(i));
                theoreticalPeakLeft = theoreticalPeakMass - MassError.getDaltonError(Properties.fragmentTolerance, theoreticalPeakMass);
                theoreticalPeakRight = theoreticalPeakMass + MassError.getDaltonError(Properties.fragmentTolerance, theoreticalPeakMass);
                maxIntensity = -1;
                for (Peak peak : peaks) {
                    peakMass = peak.getMass();
                    peakIntensity = peak.getIntensity();
                    if (peakMass >= theoreticalPeakLeft && peakMass <= theoreticalPeakRight) {
                        if (maxIntensity < peakIntensity) {
                            maxIntensity = peakIntensity;
                            peakMatchWithGreatestIntensity = peakMass;
                        }
                    }
                }
                if (maxIntensity > 0) {
                    y_ionCount++;
                    y_ionDifferenceSum += (peakMatchWithGreatestIntensity - theoreticalPeakMass);
                }
            }

            /* b-ion  */
            theoreticalPeakMass = Properties.leftIonDifference;
            for (i = 0; i < peptideLengthMinusOne; i++) {
                theoreticalPeakMass += AminoAcids.getWeightMono(peptideString.charAt(i));
                theoreticalPeakLeft = theoreticalPeakMass - MassError.getDaltonError(Properties.fragmentTolerance, theoreticalPeakMass);
                theoreticalPeakRight = theoreticalPeakMass + MassError.getDaltonError(Properties.fragmentTolerance, theoreticalPeakMass);
                maxIntensity = -1;
                for (Peak peak : peaks) {
                    peakMass = peak.getMass();
                    peakIntensity = peak.getIntensity();
                    if (peakMass >= theoreticalPeakLeft && peakMass <= theoreticalPeakRight) {
                        if (maxIntensity < peakIntensity) {
                            maxIntensity = peakIntensity;
                            peakMatchWithGreatestIntensity = peakMass;
                        }
                    }
                }
                if (maxIntensity > 0) {
                    b_ionCount++;
                    b_ionDifferenceSum += (peakMatchWithGreatestIntensity - theoreticalPeakMass);
                }
            }

        }

        U.p("Properties.rightIonDifference += " + (y_ionDifferenceSum / y_ionCount) + ";");
        U.p("Properties.leftIonDifference += " + (b_ionDifferenceSum / b_ionCount) + ";");
    }


    public void printPrecursorDifferenceToPredictedMass() {
        Match.setSortParameter(Match.SORT_BY_SPECTRUM_PEPTIDE_MASS_DIFFERENCE);
        Collections.sort(correctMatches);
        for (Match match : correctMatches) {
            double difference = match.getSpectrum().getMass() - match.getPeptide().getMass();
            U.p(match.getSpectrum().getFile().getName() + ", " + match.getSpectrum().getMass() + ": " + difference);
        }
    }

    public int countBadPrecursors() {
        int total = 0;
        for (Match match : correctMatches) {
            double difference = match.getSpectrum().getMass() - match.getPeptide().getMass();
            if (Math.abs(difference) > 2.0) total++;
        }
        return total;
    }

    public String getTestName() {
        return testName;
    }

    private ArrayList<Match> loadCorrectMatches() {
        ArrayList<Match> correctMatches = new ArrayList<Match>();
        for (Spectrum spectrum : spectra) {
            //find the file for the correct peptide
            File spectrumFile = spectrum.getFile();
            File testFolder = spectrumFile.getParentFile().getParentFile();
            File peptideFolder = new File(testFolder, "peptides");
            File peptideFile = new File(peptideFolder, spectrumFile.getName());

            //load in the correct peptide string
            boolean validPeptideFile = true;
            String correctAcidSequence = "";
            try {
                BufferedReader br = new BufferedReader(new FileReader(peptideFile));
                //read the first line;
                correctAcidSequence = br.readLine();
                //close;
                br.close();
            } catch (FileNotFoundException e) {
                validPeptideFile = false;
                e.printStackTrace();
            } catch (IOException e) {
                validPeptideFile = false;
                e.printStackTrace();
            }

            //testing that we've got a valid peptide file
            if (correctAcidSequence == null) {
                validPeptideFile = false;
            }
            correctAcidSequence = correctAcidSequence.trim();
            if (correctAcidSequence.equals("")) {
                validPeptideFile = false;
            }
            //adding to the array list
            if (validPeptideFile) {
                correctMatches.add(Properties.matchConstructor.createMatch(new MatchesSpectrum(spectrum), new Peptide(correctAcidSequence)));
            }

        }
        return correctMatches;
    }

}
