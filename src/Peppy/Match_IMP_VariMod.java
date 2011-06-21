package Peppy;

import Utilities.U;

public class Match_IMP_VariMod extends Match {
	
	//Peptide/Spectrum difference in mass
	double modificationMass;
	int modificationIndex = 0;
	
	//this is which amino acid we think has the modification
	boolean [] potentialPlacesForModification;
	
	public Match_IMP_VariMod(Spectrum spectrum, Peptide peptide) {
		setSpectrum(spectrum);
		setPeptide(peptide);
		calculateScore();
	}
	
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
		Match_IMP_VariMod match = new Match_IMP_VariMod(spectrum, peptide);
		U.p("raw score: " + match.getScore());
		
//		try {
//			SpectralVisualizerPTM.drawDeluxSpectrum(spectrum, peptide, new File ("spectrumPTM.jpg"), difference, 9);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
		
		//scores when we iterate where the mod is occuring
		double imp;
		for (int i= 0; i < acidString.length(); i++) {
			match = new Match_IMP_VariMod(spectrum, peptide);
			imp = match.calculateIMP(match.modificationMass, i);
			U.p(i + " " + acidString.charAt(i) + ": " + imp);
		}
	}
	

	
	public void calculateScore() {
		modificationMass = spectrum.getMass() - peptide.getMass();
		score = -Math.log(calculateIMP());
	}
	
	public String getScoringMethodName() {return "IMP VariMod";}
	
	public double calculateIMP() {
		if (impValue < 0) {
			double impValue1 = calculateIMP(modificationMass, 0);
			double impValue2 = calculateIMP(modificationMass, peptide.getAcidSequence().length - 1);
			if (impValue1 < 1.0E-9 || impValue2 < 1.0E-9) {
				double bestIMP = impValue1;
				double tempIMP;
				for (int i = 1; i < peptide.getAcidSequence().length - 1; i++) {
					tempIMP = calculateIMP(modificationMass, i);
					if (tempIMP < bestIMP) {
						bestIMP = tempIMP;
						modificationIndex = i;
					}
				}
				if (impValue2 < bestIMP) {
					bestIMP = impValue2;
					modificationIndex =  peptide.getAcidSequence().length - 1;
				}
				impValue = bestIMP;
			} else {
				impValue = impValue1;
				if (impValue2 < impValue) impValue = impValue2;
			}
			
		}
		return impValue;
	}
	
	
	/**
	 * This eventually calls Match's calculateIMP, but first calculates b and y ions
	 * given the modification.
	 * 
	 * @param offset
	 * @param modifiedIndex
	 * @return
	 */
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
		
		
		return calculateIMP(peptideLengthMinusOne, yIonMatchesWithHighestIntensity, bIonMatchesWithHighestIntensity);
	}
	
	public double getModificationMass() {return modificationMass;}
	
	public int getModificationIndex() {
		return modificationIndex;
	}

	public boolean hasModification() {
		return true;
	}

}
