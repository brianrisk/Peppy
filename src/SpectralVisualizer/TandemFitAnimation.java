package SpectralVisualizer;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import Peppy.AminoAcids;
import Peppy.Definitions;
import Peppy.Peak;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.Spectrum;

public class TandemFitAnimation {
	
	private static final Color yIonColor = Color.red;
	private static final  Color bIonColor = Color.blue;
	private static final  Color noIonColor = Color.gray;
	private static final  Color bothIonColor = Color.green;
	
	private  Spectrum spectrum;
	private  Peptide peptide;
	private  File folder;
	private static int fileIndex = 0;
	
	public static void main(String arenwln[]) {
		new TandemFitAnimation(
				"/Users/risk2/Documents/workspace/JavaGFS/spectra/ecoli579/E16_MSMS_1561.9093_10.t2d.txt",
				"DATGIDPVSLIAFDK"
		);
	}

	public TandemFitAnimation(String spectrumName, String peptideSequence) {
		folder = new File("animation");
		folder.mkdirs();

		spectrum = new Spectrum(spectrumName);
		peptide = new Peptide(peptideSequence);
		
		//create animation
		markMatchingIons();
		
	}
	
	private void hilightPeakAtIndex(int index) {
		ArrayList<Peak> peaks = spectrum.getPeaks();
		for (Peak peak: peaks) {
			peak.setHilighted(false);
		}
		peaks.get(index).setHilighted(true);
	}
	
	private void drawSpectrum(String message1, String message2) {
		File file = new File(folder, fileIndex + ".jpg");
		try {
			SpectralVisualizer.drawSpectrum(spectrum, 1000, 200, file, false, message1, message2);
		} catch (IOException e) {
			e.printStackTrace();
		}
		fileIndex++;
	}
	
	public void markMatchingIons() {
		
		byte [] acidSequence = peptide.getAcidSequence();
		String acidString = AminoAcids.getStringForByteArray(acidSequence);
		
		//initialize all peaks to be gray
		ArrayList<Peak> peaks = spectrum.getPeaks();
		for (Peak peak: peaks) {
			peak.setColor(noIonColor);
		}

		int i;
		double theoreticalPeakMass, peakMass;
		int peakIndex, seqIndex;
		
		//we want -1 because most of these spectra will have a match with 
		//the last theoretical peak
		int peptideLengthMinusOne = acidSequence.length - 1;
		
		double [] bIonMatchesWithHighestIntensity = new double[acidSequence.length];
		for (i = 0; i < acidSequence.length; i++) bIonMatchesWithHighestIntensity[i] = 0.0;
		double [] yIonMatchesWithHighestIntensity = new double[acidSequence.length];
		for (i = 0; i < acidSequence.length; i++) yIonMatchesWithHighestIntensity[i] = 0.0;

		
		//find the ranges around our theoretical peptides where we
		//count spectrum peaks
		double [] theoreticalPeaksLeft = new double[peptideLengthMinusOne];
		double [] theoreticalPeaksRight = new double[peptideLengthMinusOne];
		
		
		/* y-ion  */
		//computing the left and right boundaries for the ranges where our peaks should land
		theoreticalPeakMass = peptide.getMass() + Properties.rightIonDifference;
		for (i = 0; i < peptideLengthMinusOne; i++) {
			theoreticalPeakMass -= AminoAcids.getWeightMono(acidSequence[i]);
			theoreticalPeaksLeft[i] = theoreticalPeakMass - Properties.peakDifferenceThreshold;
			theoreticalPeaksRight[i] = theoreticalPeakMass + Properties.peakDifferenceThreshold;
		}
		
		peakIndex = spectrum.getPeakCount() - 1;
		seqIndex = 0;
		while (peakIndex >= 0) {
			hilightPeakAtIndex(peakIndex);
			drawSpectrum(acidString, "Y: " + acidString.substring(seqIndex + 1));
			peakMass = spectrum.getPeak(peakIndex).getMass();
			while (peakMass < theoreticalPeaksLeft[seqIndex]) {
				seqIndex++;
				if (seqIndex == peptideLengthMinusOne) break;
			}
			if (seqIndex == peptideLengthMinusOne) break;
			if (peakMass >= theoreticalPeaksLeft[seqIndex] && peakMass <= theoreticalPeaksRight[seqIndex]) {
				if (yIonMatchesWithHighestIntensity[seqIndex] < spectrum.getPeak(peakIndex).getIntensity()) {
					spectrum.getPeak(peakIndex).setColor(yIonColor);
					drawSpectrum(acidString, "Y: " + acidString.substring(seqIndex + 1));
				}
			}
			
			peakIndex--;
		}
			
		
		/* b-ion  */
		theoreticalPeakMass = Properties.leftIonDifference;
		for (i = 0; i < peptideLengthMinusOne; i++) {
			theoreticalPeakMass += AminoAcids.getWeightMono(acidSequence[i]);
			theoreticalPeaksLeft[i] = theoreticalPeakMass - Properties.peakDifferenceThreshold;
			theoreticalPeaksRight[i] = theoreticalPeakMass + Properties.peakDifferenceThreshold;
		}
		
		peakIndex = 0;
		seqIndex = 0;
		while (peakIndex < spectrum.getPeakCount()) {
			hilightPeakAtIndex(peakIndex);
			drawSpectrum(acidString, "B: " + acidString.substring(0, seqIndex + 1));
			peakMass = spectrum.getPeak(peakIndex).getMass();
			while (peakMass > theoreticalPeaksRight[seqIndex]) {
				seqIndex++;
				if (seqIndex == peptideLengthMinusOne) break;
			}
			if (seqIndex == peptideLengthMinusOne) break;
			if (peakMass >= theoreticalPeaksLeft[seqIndex] && peakMass <= theoreticalPeaksRight[seqIndex]) {
				if (bIonMatchesWithHighestIntensity[seqIndex] < spectrum.getPeak(peakIndex).getIntensity()) {
					if (spectrum.getPeak(peakIndex).getColor().equals(yIonColor)) {
						spectrum.getPeak(peakIndex).setColor(bothIonColor);						
					} else {
						spectrum.getPeak(peakIndex).setColor(bIonColor);
					}
					drawSpectrum(acidString, "B: " + acidString.substring(0, seqIndex + 1));
				}
			}
			
			peakIndex++;
		}
		
	}

}
