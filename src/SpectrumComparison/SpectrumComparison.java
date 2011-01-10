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
		Peak peak1;
		Peak peak2;
		int peaksInCommon = 0;
		double presentMassDelta;
		double presentIntensityDelta;
		// the minimum delta should never be less than zero.  It's a distance.  Duh!
		double minimumIntensityDelta = -1;
		double summedSquaredDistance = 0.0;
		while (index1 < peaks1.size() && index2 < peaks2.size()) {
			peak1 = peaks1.get(index1);
			peak2 = peaks2.get(index2);
			presentMassDelta = Math.abs(peak1.getMass() - peak2.getMass());
			if (presentMassDelta < delta ) {
				presentIntensityDelta = Math.abs(peak1.getIntensity() - peak2.getIntensity());
				//if minimumIntensityDelta has yet to be initialized
				//or it is greater than the present delta
				if (minimumIntensityDelta < 0) {
					peaksInCommon++;
					minimumIntensityDelta = presentIntensityDelta;
				} else {
					if (minimumIntensityDelta > presentIntensityDelta) {
						minimumIntensityDelta = presentIntensityDelta;
					}
				}
				/*
				 * I'm advancing the second index.  In this way I am comparing spectrum
				 * 2 to spectrum 1.  It is not the same as comparing spectrum 1 to spectrum 2.
				 * 
				 * Ideally the two would be the same.  also what would be nice is when a delta range
				 * is entered then we search for the two largest peaks between the two spectra.
				 * perhaps later we'll see how that works.
				 */
				index2++;
			} else {
				//if minimumIntensityDelta is not -1 then that means the 
				//we just now left the delta range
				if (minimumIntensityDelta >= 0) {
					summedSquaredDistance += minimumIntensityDelta * minimumIntensityDelta;
				}
				// making sure this is initialized
				minimumIntensityDelta = -1; 
				if (peak1.getMass() < peak2.getMass()) {
					index1++;
				} else {
					index2++;
				}
			}
		}
		if (peaksInCommon == 0) {
			distance = Double.MAX_VALUE;
		} else {
			distance = summedSquaredDistance / peaksInCommon;
		}
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
		double xMean = spectrum1.getAverageIntensity();
		double yMean = spectrum2.getAverageIntensity();
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
		return spectrumPeptidePair1.getPeptide().equals(spectrumPeptidePair2.getPeptide());
	}
	
	public double getDistance() {return distance;}
	
	public String toString() {
		return spectrumPeptidePair1.getPeptide().getAcidSequence() + ", " + spectrumPeptidePair2.getPeptide().getAcidSequence() + ": " + distance;
		}
	
	public Spectrum getSpectrumA() {return spectrum1;}
	public Spectrum getSpectrumB() {return spectrum2;}
	public SpectrumPeptidePair getSpectrumPeptidePairA() { return spectrumPeptidePair1;}
	public SpectrumPeptidePair getSpectrumPeptidePairB() { return spectrumPeptidePair2;}


	public int compareTo(SpectrumComparison o) {
		//sorting largest to smallest
//		if (distance > o.getDistance()) return -1;
//		if (distance < o.getDistance()) return  1;
		
		//sorting smallest to largest
		if (distance < o.getDistance()) return -1;
		if (distance > o.getDistance()) return  1;
		return 0;
	}
	

}
