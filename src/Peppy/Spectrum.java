package Peppy;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Collections;

import Math.EValueCalculator;
import Math.HasValue;
import Math.MassError;

/**
 * Holds peaks and meta-data for a MS/MS peak file.
 * Has many methods for importing peak files.  Has some analysis
 * methods like finding the average peak intensity, etc.
 * @author Brian Risk
 *
 */
public class Spectrum implements Comparable<Spectrum>, HasValue {

	private ArrayList<Peak> peaks;
	private double maxMass;
	private double mass;
	private double precursorMZ;
	private int id;
	private int charge = 0;
	private File file;
	private String MD5 = null;
	
	//E-Values
	private EValueCalculator eValueCalculator = new EValueCalculator();
	
	//in reality these values should never be less than a positive number
	private double averageIntensity = -1; 
	private double maxIntensity = -1; 
	private double medianIntensity = -1; 
	private double intensity25Percent = -1;
	private double intensity12Percent = -1;
	private double intensity06Percent = -1;
	private double minimumIntensity = -1; 
	private double coverage = -1;
	
	
	private static int sortTracker = 0;
	public final static int SORT_BY_MASS = sortTracker++;
	public final static int SORT_BY_ID = sortTracker++;
	
	//default is that we sort matches by MASS
	private static int sortParameter = SORT_BY_MASS;
	
	public Spectrum() {
		peaks = new ArrayList<Peak>();
	}
	
	public Spectrum (String s) {
		this(new File(s));
	}
	
	public Spectrum(File file) {
		this();
		this.file = file;
		Spectrum spectrum = loadSpectra(file).get(0); 
		peaks = spectrum.getPeaks();
		mass = spectrum.getMass();
	}
	
	public Spectrum(String line, File file) {
		this();
		this.file = file;
		addPrecursorFromString(line);
	}
	
	private void addPeakFromString(String s) {
		try {
			Peak p = new Peak(s);
			if (p.getMass() <= mass) peaks.add(new Peak(s));
		} catch (Exception e) {
			//don't add it if it's bad!
		}	
	}
	
	private void addPrecursorFromString(String s) {
		String [] chunks;
		chunks = s.split("\\s+"); //split on all white space
		try {
		precursorMZ = Double.parseDouble(chunks[0]);
		} catch (NumberFormatException nfe) {
			U.p("bad file: " + file.getPath());
			System.exit(1);
		}
		mass = precursorMZ - Definitions.HYDROGEN_MONO;
		if (file.getName().endsWith(".dta")) {
			charge = Integer.parseInt(chunks[1]);
		}
		//assumes txt files are really pkl
		if (file.getName().endsWith(".pkl") || file.getName().endsWith(".txt")) {
			charge = Integer.parseInt(chunks[2]);
			mass *= charge;
		}
	}
	
	/* returns a double, but we're really interested in setting the max */
	public double getMaxMass() {
		if (maxMass < 0) {
			for (Peak peak: peaks) {
				if (peak.getMass() > maxMass) maxMass = peak.getMass();
			}
		}
		return maxMass;
	}
	
	
	private void sortPeaksByIntensity() {
		Peak p;
		for (int i = 0; i < peaks.size(); i++) {
			p = (Peak) peaks.get(i);
			p.setCompareByIntensity();
		}
		Collections.sort(peaks);
	}
	
	public void sortPeaksByMass() {
		Peak p;
		for (int i = 0; i < peaks.size(); i++) {
			p = (Peak) peaks.get(i);
			p.setCompareByMass();
		}
		Collections.sort(peaks);
	}

	
	/**
	 * If the average has not been calculated already, it calculates it
	 * @return
	 */
	public double getAverageIntensity() {
		if (averageIntensity < 0) {
			averageIntensity = 0.0;
			for (Peak peak: peaks) {
				averageIntensity += peak.getIntensity();
			}
			averageIntensity /= peaks.size();
		}
		return averageIntensity;
	}
	
	/**
	 * If the maximum has not been found already, it finds it
	 * @return
	 */
	public double getMaxIntensity() {
		if (maxIntensity < 0) {
			for (Peak peak: peaks) {
				if (peak.getIntensity() > maxIntensity) maxIntensity = peak.getIntensity();
			}
		}
		return maxIntensity;
	}
	
	/**
	 * If the minimum has not been found already, it finds it
	 * @return
	 */
	public double getMinimumIntensity() {
		if (minimumIntensity < 0) {
			double min = Double.MAX_VALUE;
			for (Peak peak: peaks) {
				if (peak.getIntensity() < min) min = peak.getIntensity();
			}
			minimumIntensity = min;
		}
		return minimumIntensity;
	}
	
	/**
	 * If the median has not been found already, it finds it
	 * @return
	 */
	public double getMedianIntensity() {
		if (medianIntensity < 0) calculateDistributions();
		return medianIntensity;
	}
	
	public double getIntensity25Percent() {
		if (intensity25Percent < 0) calculateDistributions();
		return intensity25Percent;
	}
	
	public double getIntensity12Percent() {
		if (intensity12Percent < 0) calculateDistributions();
		return intensity12Percent;
	}
	
	public double getIntensity06Percent() {
		if (intensity06Percent < 0) calculateDistributions();
		return intensity06Percent;
	}
	
	private void calculateDistributions(){
		sortPeaksByIntensity();
		medianIntensity = getPeak(peaks.size() / 2).getIntensity();
		intensity25Percent = getPeak((int) (peaks.size() * .75)).getIntensity();
		intensity12Percent = getPeak((int) (peaks.size() * .875)).getIntensity();
		intensity06Percent = getPeak((int) (peaks.size() * .9375)).getIntensity();
		sortPeaksByMass();
	}
	
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
			coverage = covered / mass;
		}
		return coverage;
	}

	public Peak getPeak(int i) {return (Peak) peaks.get(i);}
	
	public ArrayList<Peak> getPeaks() {return peaks;}
	
	public double getMass() {return mass;}
	public double getPrecursorMZ() {return precursorMZ;}
	
	public int getPeakCount() {return peaks.size();}


	public int getCharge() {
		return charge;
	}

	public int getId() {
		return id;
	}

	public void setId(int id) {
		this.id = id;
	}

	public File getFile() {
		return file;
	}

	public EValueCalculator getEValueCalculator() {
		return eValueCalculator;
	}

	private void cleanPeaks() {
		//before we mess with the peak data, let's make sure we have the MD5
		MD5 = getMD5();
		
		if (Properties.highIntensityCleaning) {
			cleanPeaksKeepingHighIntensity();
		} else {
			//cleanWithWindow();
			keepStrongestPeakInRegions();
		}
		sortPeaksByMass();
		
//		markTwinPeaks();
//		markPeaksWithHighestIntensity();
		
//		keepOnlyUsedPeaks();
		
		//take exponent of peak intensities
//		if (Properties.defaultScore == Properties.DEFAULT_SCORE_TANDEM_FIT) {
//			for (Peak peak: peaks) {
//				peak.setIntensity(Math.pow(peak.getIntensity(), Properties.peakIntensityExponent));
//			}
//		}
	}
	
	/**
	 * This method is employed after other method mark
	 * peaks as "used".
	 */
	private void keepOnlyUsedPeaks() {
		for (int i = 0; i < peaks.size(); i++) {
			if (!peaks.get(i).used) {
				peaks.remove(i);
				i--;
			}
		}
	}
	
	/**
	 * The owls are not what they seem.
	 * Or, in this case, we pretend every peak is
	 * a b ion and look at all other peaks for
	 * a probable y ion.  If they exist, then both
	 * are kept.
	 */
	private void markTwinPeaks() {
		int stop = peaks.size();
		double massLowerBound, massUpperBound, perfectY;
		boolean found;
		Peak peakB, peakY;
		for (int b = 0; b < stop - 1; b++) {
			peakB = peaks.get(b);
			if (peakB.used) continue;
			
			//found starts as false because we haven't found a match yet
			found = false;
			
			//Determining the y ion
			//Start with the precursor mass
			perfectY = getMass();
			//subtract the mass of the b ion
			perfectY -= peakB.getMass();
			//since b ion is the sum of the amino acids + left ion difference,
			//that means we just subtracted the left ion difference. add that back.
			//this should now be the raw summation of the amino acids of the y ion
			perfectY += Properties.leftIonDifference;
			//add the right ion difference to get the y ion mass
			perfectY += Properties.rightIonDifference;
			massLowerBound = perfectY - MassError.getDaltonError(Properties.fragmentTolerance, perfectY);
			massUpperBound = perfectY + MassError.getDaltonError(Properties.fragmentTolerance, perfectY);
			for (int y = b + 1; y < stop; y++) {
				peakY = peaks.get(y);
				if (peakY.used) continue;
				if (peakY.getMass() > massLowerBound && peakY.getMass() < massUpperBound) {
					peakY.used = true;
					found = true;
					break;
				}
			}
			
			//handle if we found a Y match
			if (found == true) {
				peakB.used = true;
				continue;
			}
		}
	}
	
	private void keepStrongestPeakInRegions() {
		sortPeaksByIntensity();
		int start = peaks.size() - 1;
		int stop = 0;
		double lowerBound, upperBound;
		
		//note: > stop as we save final one for inside loop
		for (int i = start; i > stop; i--) {
			lowerBound = peaks.get(i).getMass() - MassError.getDaltonError(Properties.fragmentTolerance, peaks.get(i).getMass());
			upperBound = peaks.get(i).getMass() + MassError.getDaltonError(Properties.fragmentTolerance, peaks.get(i).getMass());
			for (int j = i - 1; j >= stop; j--) {
				if (peaks.get(j).getMass() > lowerBound) {
					if (peaks.get(j).getMass() < upperBound) {
						peaks.remove(j);
						i--;
						j--;
					}
				}
			}
		}
	}
	

	
	//mark the 100 most intense peaks.  Intense!
	private void cleanPeaksKeepingHighIntensity() {
		if (peaks.size() <= Properties.numberOfHighIntensityPeaksToRetain) return;
		sortPeaksByIntensity();
		ArrayList<Peak> retPeaks = new ArrayList<Peak>();
		for (int i = peaks.size() - 1; i >= peaks.size() - Properties.numberOfHighIntensityPeaksToRetain; i--) {
			retPeaks.add(peaks.get(i));
		}
		
		peaks = retPeaks;	
	}
	

	
	public static ArrayList<Spectrum> loadSpectra() {
		ArrayList<Spectrum> spectra = new ArrayList<Spectrum>();
		if (Properties.spectraDirectoryOrFile.isFile()) {
			spectra =  loadSpectra(Properties.spectraDirectoryOrFile);
		} else {
			loadSpectraFromFolder(Properties.spectraDirectoryOrFile, spectra );
		}
		
		/* remove spectra with small amount of peaks */
		int removeTally = 0;
		for (int i = 0; i < spectra.size(); i++) {
			if (spectra.get(i).getPeakCount() < Properties.minimumNumberOfPeaksForAValidSpectrum) {
				spectra.remove(i);
				i--;
				removeTally++;
			}
		}
		if (removeTally > 0) {
			U.p("removed " + removeTally + " spectra with less than " + Properties.minimumNumberOfPeaksForAValidSpectrum + " peaks");
		}
		
		/* set spectra id */
		for (int i = 0; i < spectra.size(); i++) {
			spectra.get(i).setId(i);
		}
		return spectra;
	}
	
	/**
	 * a utility method to load in a DTA/PK file
	 * @param fileName
	 * @return A ArrayList of all the spectra
	 */
	public static ArrayList<Spectrum> loadSpectra(String fileName) {
		File inFile = new File(fileName);
		return loadSpectra(inFile);
	}
	
	/**
	 * a utility method to load in a DTA/PKL file.  There may be more than
	 * one spectrum in a file, so this method returns and ArrayList
	 * @param inFile
	 * @return A ArrayList of all the spectra
	 */
	public static ArrayList<Spectrum> loadSpectra(File inFile) {
		ArrayList<Spectrum> spectra = new ArrayList<Spectrum>();
		try {
			BufferedReader inBR = new BufferedReader(new FileReader(inFile));
			String line = inBR.readLine();
			boolean fileOpen = false;
			Spectrum spectrum = null;
			while (line != null) {
				if (line.trim().equals("")) {
					if (fileOpen) {
						fileOpen = false;
						/* clean peaks also sorts the peaks by mass */
						spectrum.cleanPeaks();
						if (Properties.ignoreSpectraWithChargeGreaterThanTwo && spectrum.getCharge() <= 2) {
							spectra.add(spectrum);
						} else {
							spectra.add(spectrum);
						}
					}
				} else {
					if (fileOpen) {
						spectrum.addPeakFromString(line);
					} 
					//else this is the first line of a new spectrum
					else {
						fileOpen = true;
						spectrum = new Spectrum(line, inFile);
					}

				}
				line = inBR.readLine();
			}
			if (fileOpen) {
				/* clean peaks also sorts the peaks by mass */
				spectrum.cleanPeaks();
				if (Properties.ignoreSpectraWithChargeGreaterThanTwo && spectrum.getCharge() <= 2) {
					spectra.add(spectrum);
				} else {
					spectra.add(spectrum);
				}
			}
			inBR.close();
		}
		catch (IOException fnfe) {fnfe.printStackTrace(); U.p(inFile.getAbsolutePath());}
		catch (Exception e) {U.p(inFile.getName()); e.printStackTrace(); System.exit(1);}
		return spectra;
	}
	

	
	/**
	 * recursively goes through folder.
	 * finds all files that end in .dta or .pkl
	 * extracts the spectra from those files
	 * adds all those spectra to one big ArrayList and returns that.
	 * @param folder
	 * @param spectra
	 * @return
	 */
	public static void loadSpectraFromFolder(File folder, ArrayList<Spectrum> spectra) { 
		File [] files = folder.listFiles();
		for (int i = 0; i < files.length; i++) {
			if (files[i].isHidden()) continue;
			if (files[i].isDirectory()) {
				loadSpectraFromFolder(files[i], spectra);
				continue;
			}
			String fileName = files[i].getName().toLowerCase();
			if (fileName.endsWith(".dta") || fileName.endsWith(".pkl") || fileName.endsWith(".txt")) {
				spectra.addAll(loadSpectra(files[i]));
			}
		}
		//giving all spectra an ID
		for (int i = 0; i < spectra.size(); i++) {
			spectra.get(i).setId(i);
		}
	}
	
	/**
	 * This is slightly more general than the method it overloads.  It creates 
	 * a spectrum ArrayList which it passes to the other loadSpectraFromFolder
	 * @param fileName
	 * @return
	 */
	public static ArrayList<Spectrum> loadSpectraFromFolder(String fileName) {
		File inFile = new File(fileName);
		return loadSpectraFromFolder(inFile);
	}
	
	public static ArrayList<Spectrum> loadSpectraFromFolder(File inFile) {
		ArrayList<Spectrum> spectra = new ArrayList<Spectrum>();
		loadSpectraFromFolder(inFile, spectra);
		return spectra;
	}
	

	/**
	 * recursively goes through folder.
	 * finds all files that end in .dta or .pkl
	 * returns an array list of all of those files
	 * @param folder
	 * @param spectraFiles
	 */
	public static void loadSpectraFilesFromFolder(File folder, ArrayList<File> spectraFiles) { 
		File [] files = folder.listFiles();
		for (int i = 0; i < files.length; i++) {
			if (files[i].isHidden()) continue;
			if (files[i].isDirectory()) {
				loadSpectraFilesFromFolder(files[i], spectraFiles);
				continue;
			}
			String fileName = files[i].getName().toLowerCase();
			if (fileName.endsWith(".dta") || fileName.endsWith(".pkl") || fileName.endsWith(".txt")) {
				spectraFiles.add(files[i]);
			}
		}
	}
	

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
	}
	
	public static void setSortParameter(int sortParameter) {
		Spectrum.sortParameter = sortParameter;
	}

	public void clearEValues() {
		eValueCalculator = new EValueCalculator();
	}

	
	public double getEValue(double score) {
		return eValueCalculator.calculateEValueOfScore(score);
	}
	
	/* Assuming: https://github.com/giddingslab/peppy/issues/6
	 * P value is, roughly, the E value of a match 
	 * divided by the total number of peptides 
	 *  to which a spectrum was compared. */
	public double getPValue(double score) {
		return this.getEValue(score) / eValueCalculator.getNumberOfMatches();
	}
	
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
	}
	
	/**
	 * EXTREME CAUTION
	 * You'd better have a great reason to change this because
	 * even the slightest change will totally alter all future
	 * MD5 output
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
	}

	public double getValue() {
		return getMass();
	}
	
	

}

