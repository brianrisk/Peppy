package Peppy;

import java.io.File;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;

import Math.HasValue;
import Math.MassError;

/**
 * Holds peaks and meta-data for a MS/MS peak file.
 * Has many methods for importing peak files.  Has some analysis
 * methods like finding the average peak intensity, etc.
 * 
 * Rework to add in additional features to use the jmzml library to upload .mzml files 
 * 
 * Copyright 2013, Brian Risk
 * 
 * @author Brian Risk
 *
 */
public class Spectrum implements Comparable<Spectrum>, HasValue {

	private ArrayList<Peak> peaks;
	private double maxMass;
	private double mass = -1;
	private double precursorMZ = -1;
	private int id;
	private int fileLocus;
	private int charge = 1;
	private File file;
	private String title = "";
	//For more information about MD5 cryptographic hash function visit: http://en.wikipedia.org/wiki/MD5  
	private String MD5 = null;
	private boolean isValid = true;
	
	//in reality these values should never be less than a positive number
	private double averageIntensity = -1; 
	private double maxIntensity = -1; 
	private double medianIntensity = -1; 
	private double intensity25Percent = -1;
	private double intensity12Percent = -1;
	private double intensity06Percent = -1;
	private double minimumIntensity = -1; 
	private double coverage = -1;
	private double totalIntensity = -1;
	
	private int scanCount = -1;
	private double retentTime = -1;
	private double scanStartTime = -1;
	private double scanStopTime = -1;
	
	
	private static int sortTracker = 0;
	public final static int SORT_BY_MASS = sortTracker++;
	public final static int SORT_BY_ID = sortTracker++;
	
	//default is that we sort matches by MASS
	private static int sortParameter = SORT_BY_MASS;
	
	
/*Constructors*/
	/**
	 * Empty Constructor.  Just initializes variables.
	 */
	public Spectrum() {
		peaks = new ArrayList<Peak>();
	}

	/**
	 *
	 * @param s Location of file on HDD
	 */
	public Spectrum (String s) {
		this(new File(s));
	}//constructor
	
	/**
	 * 
	 * @param file
	 */
	public Spectrum(File file) {
		this();
		this.file = file;
		Spectrum spectrum = SpectrumLoader.loadSpectra(file).get(0); 
		peaks = spectrum.getPeaks();
		mass = spectrum.getMass();
	}//constructor


	/**
	 * calculateDistributions is a method that calculates statistics about the peaks of this Spectrum.
	 * 
	 * This method calculates: medianIntensity, 25 12 06 percent intensity
	 */
	public void calculateDistributions(){
		SpectrumLoader.sortPeaksByIntensity(this);
		medianIntensity = getPeak((int) (peaks.size() * 0.5)).getIntensity();
		intensity25Percent = getPeak((int) (peaks.size() * .75)).getIntensity();
		intensity12Percent = getPeak((int) (peaks.size() * .875)).getIntensity();
		intensity06Percent = getPeak((int) (peaks.size() * .9375)).getIntensity();
		SpectrumLoader.sortPeaksByMass(this);
		
		totalIntensity = 0;
		for (Peak peak: peaks) {
			totalIntensity += peak.getIntensity();
		}
	}//calculateDistributions
	
	
	

	
/*End Private methods for utility methods inside this file*/


	
	/**
	 * EXTREME CAUTION
	 * You'd better have a great reason to change this because
	 * even the slightest change will totally alter all future
	 * MD5 output
	 * 
	 * For more information about MD5 cryptographic hash function visit: http://en.wikipedia.org/wiki/MD5
	 * @return
	 */
	public String toStringForMD5() {
		StringBuffer out = new StringBuffer();
		for (Peak peak: peaks) {
			out.append(peak.getMass());
			out.append("\t");
			out.append(peak.getIntensity());
			out.append("\r");
		}
		return out.toString();
	}//toStringForMD5

	
	
	/**
	 * for some operations, setting a peak to "used" is handy.  This method
	 * quickly sets them all back so that used is false
	 */
	public void setAllPeaksToUnused() {
		for (Peak peak: peaks) {
			peak.used = false;
		}
	}
	
	
	
	
	
/*Getters and Setters*/
	
	public double getValue() {return getMass();}
	public Peak getPeak(int i) {return (Peak) peaks.get(i);}
	public ArrayList<Peak> getPeaks() {return peaks;}
	public double getMass() {return mass;}
	public double getPrecursorMZ() {return precursorMZ;}
	public int getPeakCount() {return peaks.size();}
	public int getCharge() {return charge;}
	public File getFile() {return file;}
	public int getId() {return id;}
	public int getFileLocus() {return fileLocus;}
	public void setMD5(String md5){MD5 = md5;};
	public void setPeaks(ArrayList<Peak> peaks){this.peaks = peaks;}
	public void setFile(File file){this.file = file;}
	public void setTitle(String title){this.title = title;}
	public String getTitle(){return this.title;}
	public void setPrecursorMZ(double m){precursorMZ = m;}
	public void setMass(double m){mass = m;}
	public void setCharge(int charge){ this.charge = charge;}
	public int getScanCount(){return this.scanCount;}
	public double getRetentTime(){return this.retentTime;}
	public void setScanCount(int scanCount){this.scanCount = scanCount;}
	public void setRetentTime(double retentTime){this.retentTime = retentTime;}
	public double getScanStartTime(){return this.scanStartTime;}
	public double getScanStopTime(){return this.scanStopTime;}
	public void setScanStartTime(double scanStartTime){this.scanStartTime = scanStartTime;}
	public void setScanStopTime(double scanStopTime){this.scanStopTime = scanStopTime;}
	/**
	 * Gets the MD5 of this object.
	 * For more information about MD5 cryptographic hash function visit: http://en.wikipedia.org/wiki/MD5
	 * @return
	 */
	public String getMD5() {
		if (MD5 != null) {
			return MD5;
		} else {
			String hashtext = null;
			String spectrumString = toStringForMD5();
			MessageDigest md5;
			try {
				md5 = MessageDigest.getInstance("MD5");
				md5.reset();
				md5.update(spectrumString.getBytes());
				byte[] digest = md5.digest();
				BigInteger bigInt = new BigInteger(1,digest);
				hashtext = bigInt.toString(16);
				// Now we need to zero pad it if you actually want the full 32 chars.
				while(hashtext.length() < 32 ){
					hashtext = "0"+hashtext;
				}
			} catch (NoSuchAlgorithmException e) {
				e.printStackTrace();
			}
			MD5 = hashtext;
			return MD5;
		}
	}//getMD5
	
	/**
	 * Given a mass window around each peak, the coverage is the percent of a spectrum
	 * that is covered by all of the peaks.  Like a target with holes shot out of it, what
	 * percent of the target is missing due to the holes.
	 * @return
	 */
	public double getCoverage() {
		if (coverage < 0) {
			double covered = 0;
			double upperBound, lowerBound, lastUpperBound = 0;
			for (Peak peak: peaks) {
				if (peak.getMass() < mass) {
					upperBound = peak.getMass() + MassError.getDaltonError(Properties.fragmentTolerance, peak.getMass());
					lowerBound = peak.getMass() - MassError.getDaltonError(Properties.fragmentTolerance, peak.getMass());
					if (lowerBound < lastUpperBound) {
						covered += upperBound - lastUpperBound;
					} else {
						covered += upperBound - lowerBound;
					}
					lastUpperBound = upperBound;
				}
			}
//			coverage = covered / (peaks.get(peaks.size() - 1).getMass() - peaks.get(0).getMass());
			covered /= getMass();
			coverage = covered;
		}
		return coverage;
	}//getCovereage
	
	
	public void setId(int id) {this.id = id;}
	public void setFileLocus(int fileLocus) {this.fileLocus = fileLocus;}
	public static void setSortParameter(int sortParameter) {Spectrum.sortParameter = sortParameter;}
	public boolean isValid() {return isValid;}
	public double getMedianIntensity() {if (medianIntensity < 0) calculateDistributions();return medianIntensity;}
	public double getIntensity25Percent() {if (intensity25Percent < 0) calculateDistributions();return intensity25Percent;}
	public double getIntensity12Percent() {if (intensity12Percent < 0) calculateDistributions();return intensity12Percent;}
	public double getIntensity06Percent() { if (intensity06Percent < 0) calculateDistributions(); return intensity06Percent;}
	
	
	public double getMaxMass() {
		if (maxMass < 0) {
			for (Peak peak: peaks) {
				if (peak.getMass() > maxMass) maxMass = peak.getMass();
			}
		}
		return maxMass;
	}//getMaxMass

	public double getAverageIntensity() {
		if (averageIntensity < 0) {
			averageIntensity = 0.0;
			for (Peak peak: peaks) {
				averageIntensity += peak.getIntensity();
			}
			averageIntensity /= peaks.size();
		}
		return averageIntensity;
	}//getAverageIntensity
	
	public double getMaxIntensity() {
		if (maxIntensity < 0) {
			for (Peak peak: peaks) {
				if (peak.getIntensity() > maxIntensity) maxIntensity = peak.getIntensity();
			}
		}
		return maxIntensity;
	}//getMaxIntensity
	

	public double getMinimumIntensity() {
		if (minimumIntensity < 0) {
			double min = Double.MAX_VALUE;
			for (Peak peak: peaks) {
				if (peak.getIntensity() < min) min = peak.getIntensity();
			}
			minimumIntensity = min;
		}
		return minimumIntensity;
	}//getMinimumIntensity
	
/*End Getters and Setters*/
	
/**
	 * compareTo method for this spectrum.  It will sort by Mass or ID depending on which parameter has been set.
	 */
	public int compareTo(Spectrum other) {
		if (sortParameter == SORT_BY_MASS) {
			if (mass > other.getMass()) return  1;
			if (mass < other.getMass()) return -1;
			return  0;
		}
		if (sortParameter == SORT_BY_ID) {
			if (id > other.getId()) return  1;
			if (id < other.getId()) return -1;
			return  0;
		}
		return 0;
	}//compareTo

public double getTotalIntensity() {
	return totalIntensity;
}
	
}//spectrum

