package Experimental;

import Peppy.Definitions;
import Peppy.Peak;
import Peppy.Spectrum;
import Peppy.SpectrumLoader;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Copyright 2013, Brian Risk
 * 
 * @author Brian Risk
 *
 */
public class SpectrumIonPairProbability {
	double fragmentToleranceInDaltons = 0.1;
	Spectrum spectrum;
	ArrayList<Double> masses;
	double pValue = -1;
	double score = -1;
	
	int n = 0;
	int k = 0;
	double p = 0;
	
	public static void main(String [] args) {
		ArrayList<Spectrum> spectra;
		Spectrum spectrum;
		SpectrumSequenceProbability ssp;
		
		spectra = SpectrumLoader.loadSpectra(new File("/Users/risk2/PeppyData/UNC/spectra/CPTAC/UNC QExactive compRef/WHIM16/run1/20312/UNC_P6_R_10.55683.55683.2.dta"));
		spectrum = spectra.get(0);
		
		SpectrumIonPairProbability sipp = new SpectrumIonPairProbability(spectrum);
	}
	
	public SpectrumIonPairProbability(Spectrum spectrum) {
		this(spectrum, 0.1);
	}
	
	public SpectrumIonPairProbability(Spectrum spectrum, double fragementToleranceInDaltons) {
		this.spectrum = spectrum;
		this.fragmentToleranceInDaltons = fragmentToleranceInDaltons;
		
		
		
		/* make masses list for searching */
		ArrayList<Peak> peaks = spectrum.getPeaks();
		masses = new ArrayList<Double>(peaks.size());
		for (Peak peak: peaks) {
			masses.add(new Double(peak.getMass()));
		}
		
		/* calculate p */
		p = (masses.size() * 2 * fragmentToleranceInDaltons) / spectrum.getMass();
		
//		1887.9522705078125
		
//		U.p(masses.get(-1 * Collections.binarySearch(masses, 1887.96)));
		
		calculate();
	}
	
	private void calculate() {
		n = masses.size();
		int foundIndex;
		double counterpart;
		double foundMass;
		for (double mass: masses) {
			counterpart = getCounterpartMass(mass);
			foundIndex = Collections.binarySearch(masses, counterpart);
			if (foundIndex >=0 ) {
				k++;
			} else {
				foundIndex *= -1;
				foundIndex -= 1;
				if (foundIndex < 0) continue;
				if (foundIndex >= n) continue;
				foundMass = masses.get(foundIndex);
				if (Math.abs(foundMass - counterpart) < fragmentToleranceInDaltons) {
					k++;
					continue;
				}
				
				foundIndex += 2;
				if (foundIndex < 0) continue;
				if (foundIndex >= n) continue;
				foundMass = masses.get(foundIndex);
				if (Math.abs(foundMass - counterpart) < fragmentToleranceInDaltons) {
					k++;
					continue;
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
	
	private double getCounterpartMass(double mass) {
		//NOTE:  would these two hydrogens be for +2 charge or what? 
		return spectrum.getMass() - mass + Definitions.HYDROGEN_MONO + Definitions.HYDROGEN_MONO;
	}

	public double getFragmentToleranceInDaltons() {
		return fragmentToleranceInDaltons;
	}

	public Spectrum getSpectrum() {
		return spectrum;
	}

	public ArrayList<Double> getMasses() {
		return masses;
	}

	public double getpValue() {
		return pValue;
	}

	public double getScore() {
		return score;
	}

	public int getN() {
		return n;
	}

	public int getK() {
		return k;
	}

	public double getP() {
		return p;
	}
	
	
	
}
