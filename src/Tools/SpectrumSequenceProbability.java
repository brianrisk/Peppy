package Tools;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

import org.apache.commons.math3.distribution.NormalDistribution;

import Math.SpectrumIonPairProbability;
import Navigator.Match;
import Peppy.Peak;
import Peppy.Spectrum;
import Peppy.SpectrumLoader;
import Peppy.U;

/**
 * Given an MS/MS spectrum, what is the probability
 * that an amino acid sequence can be derived?
 * 
A	Ala		71.03711	71.0788
R	Arg		156.10111	156.1875
N	Asn		114.04293	114.1038
D	Asp		115.02694	115.0886
C	Cys		103.00919	103.1388
E	Glu		129.04259	129.1155
Q	Gln		128.05858	128.1307
G	Gly		57.02146	57.0519
H	His		137.05891	137.1411
I	Ile		113.08406	113.1594
L	Leu		113.08406	113.1594
K	Lys		128.09496	128.1741
M	Met		131.04049	131.1926
F	Phe		147.06841	147.1766
P	Pro		97.05276	97.1167
S	Ser		87.03203	87.0782
T	Thr		101.04768	101.1051
W	Trp		186.07931	186.2132
Y	Tyr		163.06333	163.1760
V	Val		99.06841	99.1326
 * 
 * 
 * @author Brian Risk
 *
 */
public class SpectrumSequenceProbability {
	double fragmentToleranceInDaltons = 0.1;
	Spectrum spectrum;
	double pValue = -1;
	double score = -1;
	int deltaMin = 57;
	int deltaMax = 186;
	
	int n = 0;
	int k = 0;
	double p = 18.0 / (129.0);
	
	
	private static Hashtable<Integer, Integer> roundedAAMasses = new Hashtable<Integer, Integer>();
	static {
		roundedAAMasses.put(71, 71);
		roundedAAMasses.put(156, 156);
		roundedAAMasses.put(114, 114);
		roundedAAMasses.put(115, 115);
		roundedAAMasses.put(103, 103);
		roundedAAMasses.put(129, 129);
		roundedAAMasses.put(128, 128); //lys / Gln average
		roundedAAMasses.put(57, 57);
		roundedAAMasses.put(137, 137);
		roundedAAMasses.put(113, 113); // Ile / leu
		roundedAAMasses.put(131, 131);
		roundedAAMasses.put(147, 147);
		roundedAAMasses.put(97, 97);
		roundedAAMasses.put(87, 87);
		roundedAAMasses.put(101, 101);
		roundedAAMasses.put(186, 186);
		roundedAAMasses.put(163, 163);
		roundedAAMasses.put(99, 99);

	}
	
	
	
	public static void main(String [] args) {
		/*
		 * loading matches
		 */
		U.p("loading matches");
		ArrayList<Match> matches = new ArrayList<Match>();
		matches.addAll(Match.loadMatches(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM16-CompRef/1 WHIM16 - gencode proteome/report.txt")));
		matches.addAll(Match.loadMatches(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM16-CompRef/4 WHIM16 - mouse proteome/report.txt")));
		matches.addAll(Match.loadMatches(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM16-CompRef/6 WHIM16 - varimod/report.txt")));
		
		/*
		 * Finding the list of spectrum files
		 */
		U.p("finding spectrum files");
		
		Hashtable<String, File> usedSpectra = new Hashtable<String, File>();
		for (Match match: matches) {
			File spectrumFile = match.getFile("FilePath");
			if (spectrumFile == null) continue;
			String spectrumMD5 = match.getString("spectrumMD5");
			usedSpectra.put(spectrumMD5, spectrumFile );
		}
		ArrayList<File> spectrumFiles = new ArrayList<File>(usedSpectra.values());
		U.p("total spectrum files: " + spectrumFiles.size());
		
		/*
		 * Loading the whole mess of spectra
		 */
		U.startStopwatch();
		U.p("loading all of the spectra from the experiment");
		ArrayList<Spectrum> allSpectra = SpectrumLoader.loadSpectra(new File("/Users/risk2/PeppyData/UNC/spectra/CPTAC/UNC QExactive compRef/WHIM16"));
		U.stopStopwatch();
		U.p("this is the full spectral set size: " + allSpectra.size());
		
		/*
		 * create list of used spectra
		 */
		U.p("create list of used spectra");
		ArrayList<Spectrum> spectra = new ArrayList<Spectrum>(spectrumFiles.size());
		for (Spectrum spectrum: allSpectra) {
			if (usedSpectra.get(spectrum.getMD5()) != null) {
				spectra.add(spectrum);
			}
		}
		U.p("found this many for used: " + spectra.size());
		
		/*
		 * loading the spectra from files
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
		 */
		spectra = allSpectra;
		
		
		/*
		 * computing probabilities for spectra
		 */
		U.p("computing probablities");
		U.startStopwatch();
		Hashtable<Integer, Integer> tallies = new Hashtable<Integer, Integer>();
		int score;
		double scoreDouble;
		Integer tallyValue;
		SpectrumIonPairProbability sipp ;
		SpectrumSequenceProbability ssp;
		for (Spectrum spectrum: spectra) {
			sipp = new SpectrumIonPairProbability(spectrum);
			ssp = new SpectrumSequenceProbability(spectrum);
			scoreDouble = sipp.getScore();
//			scoreDouble += ssp.getScore();
			score = (int) Math.round(10 * scoreDouble);
			tallyValue = tallies.get(score);
			if (tallyValue == null) {
				tallyValue = 1;
			} else {
				tallyValue = tallyValue + 1;
			}
			tallies.put(score, tallyValue);
		}
		U.stopStopwatch();
		
		/*
		 * writing results
		 */
		U.p("writing results");
		try {
			PrintWriter pw = new PrintWriter(new FileWriter("probability tallies.txt"));
			ArrayList<Integer> keys = new ArrayList<Integer>(tallies.keySet());
			Collections.sort(keys);
			for (Integer key: keys) {
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
		this(spectrum, 0.1);
	}
	
	public SpectrumSequenceProbability(Spectrum spectrum, double fragmentToleranceInDaltons) {
		this.spectrum = spectrum;
		this.fragmentToleranceInDaltons = fragmentToleranceInDaltons;
		p *= fragmentToleranceInDaltons * 2;
		calculate();
	}
	
	public double getPValue() {
		return pValue;
	}
	
	private void calculate() {
		ArrayList<Peak> peaks = spectrum.getPeaks();
		double massA, massB;
		int delta;
		Peak peakA, peakB;
		for (int indexA = 0; indexA < peaks.size() - 1; indexA++) {
			peakA = peaks.get(indexA);
			massA = peakA.getMass();
			for (int indexB = indexA + 1; indexB < peaks.size(); indexB ++) {
				peakB = peaks.get(indexB);
				massB = peakB.getMass();
				delta = (int) Math.round(massB - massA);
				if (delta < deltaMin) continue;
				if (delta > deltaMax) break;
				n++;
				if (roundedAAMasses.get(delta) != null) {
//					indexA = indexB - 1;
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

	public int getK() {
		return k;
	}

	public double getScore() {
		return score;
	}

}
