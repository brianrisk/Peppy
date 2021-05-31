package Experimental;

import Graphs.Point;
import Graphs.Scatter;
import Graphs.ScatterVisualizer;
import Navigator.MatchRow;
import Peppy.*;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

/**
 * Given an MS/MS spectrum, what is the probability
 * that an amino acid sequence can be derived?
 * <p>
 * A	Ala		71.03711	71.0788
 * R	Arg		156.10111	156.1875
 * N	Asn		114.04293	114.1038
 * D	Asp		115.02694	115.0886
 * C	Cys		103.00919	103.1388
 * E	Glu		129.04259	129.1155
 * Q	Gln		128.05858	128.1307
 * G	Gly		57.02146	57.0519
 * H	His		137.05891	137.1411
 * I	Ile		113.08406	113.1594
 * L	Leu		113.08406	113.1594
 * K	Lys		128.09496	128.1741
 * M	Met		131.04049	131.1926
 * F	Phe		147.06841	147.1766
 * P	Pro		97.05276	97.1167
 * S	Ser		87.03203	87.0782
 * T	Thr		101.04768	101.1051
 * W	Trp		186.07931	186.2132
 * Y	Tyr		163.06333	163.1760
 * V	Val		99.06841	99.1326
 *
 * @author Brian Risk
 */
public class SpectrumSequenceProbability {
    static double fragmentTolerance = 0.015;
    Spectrum spectrum;
    double pValue = -1;
    double score = -1;
    double deltaMin = 57.01; // 57.02146
    double deltaMax = 186.09; //186.07931

    int n = 0;
    int k = 0;
    double p = 18.0 / (129.0);


    private char getAcid(double mass) {
        int roundedMass = (int) Math.round(mass);
        if (roundedMass == 57) {
            if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.G)) < fragmentTolerance) return 'G';
        }
        if (roundedMass == 71) {
            if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.A)) < fragmentTolerance) return 'A';
        }
        if (roundedMass == 87) {
            if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.S)) < fragmentTolerance) return 'S';
        }
        if (roundedMass == 97) {
            if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.P)) < fragmentTolerance) return 'P';
        }
        if (roundedMass == 99) {
            if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.V)) < fragmentTolerance) return 'V';
        }
        if (roundedMass == 101) {
            if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.T)) < fragmentTolerance) return 'T';
        }
        if (roundedMass == 103) {
            if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.C)) < fragmentTolerance) return 'C';
        }
        if (roundedMass == 113) {
            if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.I)) < fragmentTolerance) return 'I';
        } /* or L */
        if (roundedMass == 114) {
            if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.N)) < fragmentTolerance) return 'N';
        }
        if (roundedMass == 115) {
            if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.D)) < fragmentTolerance) return 'D';
        }
        if (roundedMass == 128) {
            if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.Q)) < fragmentTolerance) return 'Q';
        }
        if (roundedMass == 128) {
            if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.K)) < fragmentTolerance) return 'K';
        }
        if (roundedMass == 129) {
            if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.E)) < fragmentTolerance) return 'E';
        }
        if (roundedMass == 131) {
            if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.M)) < fragmentTolerance) return 'M';
        }
        if (roundedMass == 137) {
            if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.H)) < fragmentTolerance) return 'H';
        }
        if (roundedMass == 147) {
            if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.F)) < fragmentTolerance) return 'F';
        }
        if (roundedMass == 156) {
            if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.R)) < fragmentTolerance) return 'R';
        }
        if (roundedMass == 163) {
            if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.Y)) < fragmentTolerance) return 'Y';
        }
        if (roundedMass == 186) {
            if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.W)) < fragmentTolerance) return 'W';
        }
        return '*';
    }


    public static void main(String[] args) {
        /*
         * loading matches
         */
        U.p("loading matches");
        ArrayList<MatchRow> matches = new ArrayList<MatchRow>();
//		matches.addAll(MatchRow.loadMatches(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM16-CompRef/1 WHIM16 - gencode proteome/report.txt")));
//		matches.addAll(MatchRow.loadMatches(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM16-CompRef/4 WHIM16 - mouse proteome/report.txt")));
//		matches.addAll(MatchRow.loadMatches(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM16-CompRef/6 WHIM16 - varimod/report.txt")));


//		matches.addAll(MatchRow.loadMatches(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/UNC-ENCODE-WCL/1 UNC-ENCODE-WCL - UniProt_Human_2012_03.fasta/report.txt")));
//		matches.addAll(MatchRow.loadMatches(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/UNC-CPTAC-QEXATIVE/1 UNC-CPTAC-QEXATIVE - UniProt_Human_2012_03.fasta/report.txt")));

        matches.addAll(MatchRow.loadMatches(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/WASHU-CPTAC-TTOF/1 WASHU-CPTAC-TTOF - UniProt_Human_2012_03.fasta/report.txt")));
//		matches.addAll(MatchRow.loadMatches(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/KAPP/1 spectra - UniProt_Human_2012_03.fasta/report.txt")));



        /*
         * Finding the list of spectrum files
         */
        U.p("finding spectrum files");

        Hashtable<String, MatchRow> usedSpectra = new Hashtable<String, MatchRow>();
        for (MatchRow match : matches) {
            File spectrumFile = match.getFile("FilePath");
            if (spectrumFile == null) continue;
            String spectrumMD5 = match.getString("spectrumMD5");
            usedSpectra.put(spectrumMD5, match);
        }
        U.p("total spectrum files: " + usedSpectra.size());

        /*
         * Loading the whole mess of spectra
         */
        U.startStopwatch();
        U.p("loading all of the spectra from the experiment");
        File spectrumDirectory = matches.get(0).getFile("FilePath").getParentFile();
        ArrayList<Spectrum> allSpectra = SpectrumLoader.loadSpectra(spectrumDirectory);

        U.stopStopwatch();
        U.p("this is the full spectral set size: " + allSpectra.size());

        /*
         * create list of used spectra
         */
        U.p("create list of used spectra");
        ArrayList<Spectrum> spectra = new ArrayList<Spectrum>(usedSpectra.size());
        for (Spectrum spectrum : allSpectra) {
            if (usedSpectra.get(spectrum.getMD5()) != null) {
                if (usedSpectra.get(spectrum.getMD5()).getFile("FilePath") != null) {
                    spectra.add(spectrum);
                }
            }
        }
        U.p("found this many for used: " + spectra.size());

        /*
         * loading the spectra from files
         *
         * NOTE: for the time being, this is commented as, for some reason
         * this is taking a huge amount of time.
         *
         */
//		U.p("loading spectra");
//		ArrayList<Spectrum> spectra = new ArrayList<Spectrum>(spectrumFiles.size());
//		int count = 0;
//		U.startStopwatch();
//		for (File spectrumFile: spectrumFiles) {
//			spectra.addAll(SpectrumLoader.loadSpectraFromFile(spectrumFile));
//			count++;
//			if (count % 1000 == 0) {
//				U.stopStopwatch();
//				U.startStopwatch();
//				U.p(count);
//				U.p();
//			}
//		}
//		U.stopStopwatch();

        /*
         * Use this if we want to use all spectra, not just those found
         *
         * We would do this if we wanted to see properties for the full spectral
         * body, not just those that matched at the 1% FDR
         *
         */
//		spectra = allSpectra;


        /*
         * computing probabilities for spectra
         */
        U.p("computing probablities");
        U.startStopwatch();

        /* setting up variables */
        Hashtable<Integer, Integer> tallies = new Hashtable<Integer, Integer>();
        int score;
        double scoreDouble;
        Integer tallyValue;
        SpectrumIonPairProbability sipp;
        SpectrumSequenceProbability ssp;

        /* setting up writer for scatter plot */
        try {
            PrintWriter scatterWriter = new PrintWriter(new FileWriter("spectrum to score scatter.txt"));
            Scatter scatter = new Scatter();

            /* computing score  and tallies for each spectrum */
            for (Spectrum spectrum : spectra) {
                sipp = new SpectrumIonPairProbability(spectrum);
                ssp = new SpectrumSequenceProbability(spectrum);
                scoreDouble = ssp.getScore();
//				scoreDouble += sipp.getScore();
                score = (int) Math.round(10 * scoreDouble);
                tallyValue = tallies.get(score);
                if (tallyValue == null) {
                    tallyValue = 1;
                } else {
                    tallyValue = tallyValue + 1;
                }
                tallies.put(score, tallyValue);

                /* printing result for scatter */
                scatterWriter.println(usedSpectra.get(spectrum.getMD5()).getScore() + "\t" + scoreDouble);
                scatter.addPoint(new Point(usedSpectra.get(spectrum.getMD5()).getScore(), scoreDouble));
            }

            ScatterVisualizer scatterVisualizer = new ScatterVisualizer(scatter, new File("visualized score to spectrum scatter.png"));
            ;
            scatterVisualizer.draw();

            scatterWriter.flush();
            scatterWriter.close();
        } catch (IOException e1) {
            // TODO Auto-generated catch block
            e1.printStackTrace();
        }

        U.stopStopwatch();

        /*
         * writing tally results
         */
        U.p("writing results");
        try {
            PrintWriter pw = new PrintWriter(new FileWriter("probability tallies.txt"));
            ArrayList<Integer> keys = new ArrayList<Integer>(tallies.keySet());
            Collections.sort(keys);
            for (Integer key : keys) {
                pw.println(key + "\t" + tallies.get(key));
            }
            pw.flush();
            pw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        U.p("done");


//		ArrayList<Spectrum> spectra;
//		Spectrum spectrum;
//		SpectrumSequenceProbability ssp;
//		
//		spectra = SpectrumLoader.loadSpectra(new File("/Users/risk2/PeppyData/UNC/spectra/CPTAC/UNC QExactive compRef/WHIM16/run1/20312/UNC_P6_R_10.55683.55683.2.dta"));
//		spectrum = spectra.get(0);
//		ssp = new SpectrumSequenceProbability(spectrum, 0.11);
//		if (Double.isInfinite(ssp.getScore())) U.p("is infinie");
//		U.p(ssp.getScore());
//		U.p();
//		
//		spectra = SpectrumLoader.loadSpectra(new File("/Users/risk2/PeppyData/UNC/spectra/CPTAC/UNC QExactive compRef/WHIM16/run1/20327/UNC_P6_R_18.19929.19929.2.dta"));
//		spectrum = spectra.get(0);
//		ssp = new SpectrumSequenceProbability(spectrum, 0.1);
//		U.p(-Math.log(ssp.getPValue()));
//		U.p(ssp.getK());
//		U.p();
//		


    }


    public SpectrumSequenceProbability(Spectrum spectrum) {
        this.spectrum = spectrum;
        p *= fragmentTolerance * 2;
        calculate();
    }


    private void calculate() {
        ArrayList<Peak> peaks = spectrum.getPeaks();
        double massA, massB, delta;
        Peak peakA, peakB;
        for (int indexA = 0; indexA < peaks.size() - 1; indexA++) {
            peakA = peaks.get(indexA);
            massA = peakA.getMass();
            for (int indexB = indexA + 1; indexB < peaks.size(); indexB++) {
                peakB = peaks.get(indexB);
                massB = peakB.getMass();
                delta = massB - massA;
                if (delta < deltaMin) continue;
                if (delta > deltaMax) break;
                n++;
                if (getAcid(delta) != '*') {
                    k++;
                }
            }
        }

        /* if no matches at all */
        if (n == 0) {
            score = 0;
            return;
        }

        double median = p * n;
        if (k < median) {
            pValue = 1;
        } else {
            double kAdjusted = median - k;
//			double kAdjusted = k - median;
            NormalDistribution nd = new NormalDistribution(n * p, n * p * (1 - p));
            pValue = nd.cumulativeProbability(kAdjusted);
        }

        /* avoiding scores of "infinity" */
        if (pValue == 0) {
            score = 40;
        } else {
            score = -Math.log(pValue);
        }

    }

    public double getScore() {
        return score;
    }

}
