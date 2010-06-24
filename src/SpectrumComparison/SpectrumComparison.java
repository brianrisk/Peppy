package SpectrumComparison;

import java.util.ArrayList;
import Peppy.Peak;
import Peppy.Spectrum;


public class SpectrumComparison implements Comparable<SpectrumComparison> {
	
	private SpectrumPeptidePair spectrumPeptidePair1;
	private SpectrumPeptidePair spectrumPeptidePair2;
	private Spectrum spectrum1;
	private Spectrum spectrum2;
	private double delta;
	private double distance = 0.0;
	
	/*
	 * If we happen to know the sequences that begat the spectra
	 * (eg. we are performing a test on our comparison algorithm)
	 */
	private String sequence1 = "";
	private String sequence2 = "";
	
	/**
	 * @param spectrum1
	 * @param spectrum2
	 * @param delta - the delta is how far away a peak of one spectrum can be from another peak of a different spectrum and still be considered the same peak
	 */
	public SpectrumComparison(SpectrumPeptidePair spectrumPeptidePair1, SpectrumPeptidePair spectrumPeptidePair2, double delta) {
		this.spectrumPeptidePair1 = spectrumPeptidePair1;
		this.spectrumPeptidePair2 = spectrumPeptidePair2;
		this.spectrum1 = spectrumPeptidePair1.getSpectrum();
		this.spectrum2 = spectrumPeptidePair2.getSpectrum();
		this.delta = delta;
		computeDistance();
//		computeProductMomentCorrelationCoefficient();
	}
	
	public double computeDistance() {
		ArrayList<Peak> peaks1 = spectrum1.getPeaks();
		ArrayList<Peak> peaks2 = spectrum2.getPeaks();
		int index1 = 0; 
		int index2 = 0;
		Peak peak1 = peaks1.get(0);
		Peak peak2 = peaks2.get(0);
		int peaksInCommon = 0;
		double xMean = spectrum1.getCalculatedAverageIntensity();
		double yMean = spectrum2.getCalculatedAverageIntensity();	
		while (index1 < peaks1.size() && index2 < peaks2.size()) {
			peak1 = peaks1.get(index1);
			peak2 = peaks2.get(index2);
			
			if (peak1.getMass() - peak2.getMass() < delta) {
				if (peak1.getIntensity() > xMean && peak2.getIntensity() > yMean) peaksInCommon++;
				index1++;
				index2++;
			} else {
				if (peak1.getMass() < peak2.getMass()) {
					if (index1 < peaks1.size()) {
						index1++;
					}
				} else {
					index2++;
				}
			}
		}
		//distance = (double) peaksInCommon / (xMean * yMean * (peaks1.size() + peaks2.size()));
		distance = (double) peaksInCommon / (peaks1.size() + peaks2.size());
		return distance;
	}
	
	public double computeProductMomentCorrelationCoefficient() {
		ArrayList<Peak> peaks1 = spectrum1.getPeaks();
		ArrayList<Peak> peaks2 = spectrum2.getPeaks();
		int index1 = 0;
		int index2 = 0;
		Peak peak1 = peaks1.get(0);
		Peak peak2 = peaks2.get(0);
		double sigmaX1Y1 = 0.0;
		double sigmaX2 = 0.0;
		double sigmaY2 = 0.0;
		double xMean = spectrum1.getCalculatedAverageIntensity();
		double yMean = spectrum2.getCalculatedAverageIntensity();
		double xDifference, yDifference;
		while (index1 < peaks1.size() && index2 < peaks2.size()) {
			//making sure we're not out of bounds
			if (index1 < peaks1.size()) peak1 = peaks1.get(index1);
			if (index2 < peaks2.size()) peak2 = peaks2.get(index2);
			
			if (peak1.getMass() - peak2.getMass() < delta) {
				xDifference = peak1.getIntensity() - xMean;
				yDifference = peak2.getIntensity() - yMean;
				sigmaX1Y1 += xDifference * yDifference;
				sigmaX2 += xDifference * xDifference;
				sigmaY2 += yDifference * yDifference;
				index1++;
				index2++;
			} else {
				if (peak1.getMass() < peak2.getMass()) {
					xDifference = peak1.getIntensity() - xMean;
					sigmaX1Y1 += xDifference * yMean * -1;
					sigmaX2 += xDifference * xDifference;
					sigmaY2 += yMean * yMean;
					if (index1 < peaks1.size()) {
						index1++;
					}
				} else {
					yDifference = peak2.getIntensity() - yMean;
					sigmaX1Y1 += xMean * yDifference * -1;
					sigmaX2 += xMean * xMean;
					sigmaY2 += yDifference * yDifference;
					index2++;
				}
			}
		}
		distance = Math.abs(sigmaX1Y1 / Math.sqrt(sigmaX2 * sigmaY2));
		return distance;
	}
	
	public boolean isEqual() {
		return spectrumPeptidePair1.getPeptide().getAcidSequence().equals(spectrumPeptidePair2.getPeptide().getAcidSequence());
	}
	
	public double getDistance() {return distance;}
	
	public String toString() {return spectrumPeptidePair1.getPeptide().getAcidSequence() + ", " + spectrumPeptidePair2.getPeptide().getAcidSequence();}


	public int compareTo(SpectrumComparison o) {
		//sorting largest to smallest
		if (distance > o.getDistance()) return -1;
		if (distance < o.getDistance()) return  1;
		
		//sorting smallest to largest
//		if (distance < o.getDistance()) return -1;
//		if (distance > o.getDistance()) return  1;
		return 0;
	}
	

}
