package Peppy;

import java.io.File;
import java.io.IOException;

import Math.MathFunctions;
import SpectralVisualizer.SpectralVisualizerPTM;
import Utilities.U;

public class MatchPTM extends Match {
	
	//Peptide/Spectrum difference in mass
	double difference;
	
	//this is which amino acid we think has the modification
	boolean [] potentialPlacesForModification;
	
	public static void test() {
		//load our weird spectrum
//		Spectrum spectrum = Spectrum.loadSpectra("/Users/risk2/PeppyOverflow/tests/kapp-just-modified/spectra/JNA-and-JGP-AAP-27Aug03_HUPO.1948.1948.1.dta").get(0);
//		Spectrum spectrum = Spectrum.loadSpectra("/Users/risk2/PeppyOverflow/tests/kapp-just-modified/spectra/JNA-and-JGP-AAP-27Aug03_HUPO.1987.1987.3.dta").get(0);
//		Spectrum spectrum = Spectrum.loadSpectra("/Users/risk2/PeppyOverflow/tests/kapp-just-modified/spectra/JNA-and-JGP-AAP-27Aug03_HUPO.2616.2616.3.dta").get(0);
		Spectrum spectrum = Spectrum.loadSpectra("/Users/risk2/PeppyOverflow/tests/kapp-just-modified/spectra/JNA-and-JGP-AAP-27Aug03_HUPO.2695.2695.2.dta").get(0);
//		Spectrum spectrum = Spectrum.loadSpectra("/Users/risk2/PeppyOverflow/tests/kapp-just-modified/spectra/JNA-and-JGP-AAP-27Aug03_HUPO.2856.2856.3.dta").get(0);
		
		
		//make our peptide
//		String acidString = "AVaDDFAAFVEK";
//		String acidString = "SVSGKPQYaVLVPSLLHTETTEK";
//		String acidString = "VVSMDENFHPLNELIPLVYIQDPkGNR";
		String acidString = "EQLkAVMDDFAAFVEK";
//		String acidString = "VFDEFKPLVEEPQNLIkQNCELFEQLGEYK";
		
		
		Peptide peptide = new Peptide(acidString.toUpperCase());
		
		//quick report
		U.p("peptide mass: " + peptide.getMass());
		U.p("spectrum mass: " + spectrum.getMass());
		double difference = ( spectrum.getMass() - peptide.getMass());
		U.p("difference: " + difference);
		MatchPTM match = new MatchPTM(spectrum, peptide);
		U.p("raw score: " + match.getScore());
		
		try {
			SpectralVisualizerPTM.drawDeluxSpectrum(spectrum, peptide, new File ("spectrumPTM.jpg"), difference, 9);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//scores when we iterate where the mod is occuring
		double imp;
		for (int i= 0; i < acidString.length(); i++) {
			match = new MatchPTM(spectrum, peptide);
			imp = match.calculateIMP(match.difference, i);
			U.p(i + " " + acidString.charAt(i) + ": " + imp);
		}
	}
	
	public MatchPTM(Spectrum spectrum, Peptide peptide) {
		super(spectrum, peptide, false);
		difference = spectrum.getMass() - peptide.getMass();
		calculateScore();
	}
	
	public void calculateScore() {
		score = -Math.log(calculateIMP());
	}
	
	public double calculateIMP() {
		if (impValue < 0) {
			double impValue1 = calculateIMP(difference, 0);
			double impValue2 = calculateIMP(difference, peptide.getAcidSequence().length - 1);
			if (impValue1 < 1.0E-9 || impValue2 < 1.0E-9) {
				double bestIMP = impValue1;
				double tempIMP;
				for (int i = 1; i < peptide.getAcidSequence().length - 1; i++) {
					tempIMP = calculateIMP(difference, i);
					if (tempIMP < bestIMP) bestIMP = tempIMP; 
				}
				if (impValue2 < bestIMP) bestIMP = impValue2;
				impValue = bestIMP;
			} else {
				impValue = impValue1;
				if (impValue2 < impValue) impValue = impValue2;
			}
			
		}
		return impValue;
	}
	
	public double calculateIMP(double offset, int modifiedIndex) {
		byte [] acidSequence = peptide.getAcidSequence();

		int i;
		boolean atLeastOneMatch = false;
		double theoreticalPeakMass, peakMass;
		int peakIndex, seqIndex;
		ionMatchTally = 0;
		
		//we want -1 because most of these spectra will have a match with 
		//the last theoretical peak
		int peptideLengthMinusOne = acidSequence.length - 1;
		
		double [] bIonMatchesWithHighestIntensity = new double[peptideLengthMinusOne];
		double [] yIonMatchesWithHighestIntensity = new double[peptideLengthMinusOne];

		//find the ranges around our theoretical peptides where we
		//count spectrum peaks
		double [] theoreticalPeaksLeft = new double[peptideLengthMinusOne];
		double [] theoreticalPeaksRight = new double[peptideLengthMinusOne];
		
		
		/* y-ion  */
		//computing the left and right boundaries for the ranges where our peaks should land
		theoreticalPeakMass = offset + peptide.getMass() + Properties.rightIonDifference;
		for (i = 0; i < peptideLengthMinusOne; i++) {
			theoreticalPeakMass -= AminoAcids.getWeightMono(acidSequence[i]);
			if (i == modifiedIndex) theoreticalPeakMass -= offset;
			theoreticalPeaksLeft[i] = theoreticalPeakMass - Properties.peakDifferenceThreshold;
			theoreticalPeaksRight[i] = theoreticalPeakMass + Properties.peakDifferenceThreshold;	
		}
		
		peakIndex = spectrum.getPeakCount() - 1;
		seqIndex = 0;
		while (peakIndex >= 0) {
			peakMass = spectrum.getPeak(peakIndex).getMass();
			while (peakMass < theoreticalPeaksLeft[seqIndex]) {
				seqIndex++;
				if (seqIndex == peptideLengthMinusOne) break;
			}
			if (seqIndex == peptideLengthMinusOne) break;
			if (peakMass >= theoreticalPeaksLeft[seqIndex] && peakMass <= theoreticalPeaksRight[seqIndex]) {
				atLeastOneMatch = true;
				if (yIonMatchesWithHighestIntensity[seqIndex] < spectrum.getPeak(peakIndex).getIntensity()) {
					yIonMatchesWithHighestIntensity[seqIndex] = spectrum.getPeak(peakIndex).getIntensity();
				}
			}
			
			peakIndex--;
		}

		//if 0 matches so far, just get out.
		if (!atLeastOneMatch) {
			impValue = 1;
			return impValue;
		}
			
		
		/* b-ion  */
		theoreticalPeakMass = Properties.leftIonDifference;
		for (i = 0; i < peptideLengthMinusOne; i++) {
			theoreticalPeakMass += AminoAcids.getWeightMono(acidSequence[i]);
			if (i == modifiedIndex) theoreticalPeakMass += offset;
			theoreticalPeaksLeft[i] = theoreticalPeakMass - Properties.peakDifferenceThreshold;
			theoreticalPeaksRight[i] = theoreticalPeakMass + Properties.peakDifferenceThreshold;
		}
		
		peakIndex = 0;
		seqIndex = 0;
		while (peakIndex < spectrum.getPeakCount()) {
			peakMass = spectrum.getPeak(peakIndex).getMass();
			while (peakMass > theoreticalPeaksRight[seqIndex]) {
				seqIndex++;
				if (seqIndex == peptideLengthMinusOne) break;
			}
			if (seqIndex == peptideLengthMinusOne) break;
			if (peakMass >= theoreticalPeaksLeft[seqIndex] && peakMass <= theoreticalPeaksRight[seqIndex]) {
				atLeastOneMatch = true;
				if (bIonMatchesWithHighestIntensity[seqIndex] < spectrum.getPeak(peakIndex).getIntensity()) {
					bIonMatchesWithHighestIntensity[seqIndex] = spectrum.getPeak(peakIndex).getIntensity();
				}
			}
			
			peakIndex++;
		}
		
		//if 0 matches so far, just get out.
		if (!atLeastOneMatch) {
			impValue = 1;
			return impValue;
		}
		
		//find out ionMatchTally and intensity distributions
		int totalIonsAbove50 = 0;
		int totalIonsAbove25 = 0;
		int totalIonsAbove12 = 0;
		int totalIonsAbove06 = 0;
		double totalMatchingIntensity = 0.0;
		boolean yIonMatch, bIonMatch, ionMatch, previousIonMatch = false;
		int acidMatchTally = 0;
		int consecutiveAcidTally = 0;
		for (i = 0; i < peptideLengthMinusOne; i++) {
			yIonMatch = yIonMatchesWithHighestIntensity[i] > 0.0;
			bIonMatch = bIonMatchesWithHighestIntensity[i] > 0.0;
			ionMatch = bIonMatch || yIonMatch;
			if (ionMatch) acidMatchTally++;
			if (previousIonMatch && ionMatch) consecutiveAcidTally++;
			if (yIonMatch) {
				ionMatchTally++;
				totalMatchingIntensity += yIonMatchesWithHighestIntensity[i];
				if (yIonMatchesWithHighestIntensity[i] > spectrum.getMedianIntensity()) {
					totalIonsAbove50++;
					if (yIonMatchesWithHighestIntensity[i] > spectrum.getIntensity25Percent()) {
						totalIonsAbove25++;
						if (yIonMatchesWithHighestIntensity[i] > spectrum.getIntensity12Percent()) {
							totalIonsAbove12++;
							if (yIonMatchesWithHighestIntensity[i] > spectrum.getIntensity06Percent()) {
								totalIonsAbove06++;
							}
						}
					}
				}
			}
			if (bIonMatch) {
				ionMatchTally++;
				totalMatchingIntensity += bIonMatchesWithHighestIntensity[i];
				if (bIonMatchesWithHighestIntensity[i] > spectrum.getMedianIntensity()) {
					totalIonsAbove50++;
					if (bIonMatchesWithHighestIntensity[i] > spectrum.getIntensity25Percent()) {
						totalIonsAbove25++;
						if (bIonMatchesWithHighestIntensity[i] > spectrum.getIntensity12Percent()) {
							totalIonsAbove12++;
							if (bIonMatchesWithHighestIntensity[i] > spectrum.getIntensity06Percent()) {
								totalIonsAbove06++;
							}
						}
					}
				}
			}
			previousIonMatch = ionMatch;
			
		}
		
		
		//Variables for binomial probabilities
		int n;
		int k;
		double p;
		
		//peak match probability is the binomial distribution
		n = acidSequence.length * 2;
		k = ionMatchTally;
		p = spectrum.getCoverage();
		double peakMatchProbability = MathFunctions.getBinomialProbability(n, k, p);
		
		
		//TODO: optimize this.  it is a great place to improve performance
		if (peakMatchProbability > 0.5) {
			impValue = 1;
			return impValue;
		}
		
		//consecutive match probability
		n = peptideLengthMinusOne;
		k = consecutiveAcidTally;
		p = (double) acidMatchTally / peptideLengthMinusOne;
		double consecutiveProbability = MathFunctions.getBinomialProbability(n, k, p);
		
		
		//probability of ions being above thresholds
		double intensityProbability = 1;
		if (ionMatchTally > 0) {
			n = ionMatchTally;
			k = totalIonsAbove50;
			intensityProbability = MathFunctions.getCachedBinomialProbability50(n, k);
			if (totalIonsAbove50 > 0) {
				n = totalIonsAbove50;
				k = totalIonsAbove25;
				intensityProbability *= MathFunctions.getCachedBinomialProbability50(n, k);
				if (totalIonsAbove25 > 0) {
					n = totalIonsAbove25;
					k = totalIonsAbove12;
					intensityProbability *= MathFunctions.getCachedBinomialProbability50(n, k);
					if (totalIonsAbove12 > 0) {
						n = totalIonsAbove12;
						k = totalIonsAbove06;
						intensityProbability *= MathFunctions.getCachedBinomialProbability50(n, k);
					}
				}
			}
		}
		
		impValue =  peakMatchProbability  * intensityProbability * consecutiveProbability;
		return impValue;
	}
	
	public double getDifference() {return difference;}
	
	public boolean hasModification() {
		return true;
	}

}
