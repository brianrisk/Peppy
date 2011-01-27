package Peppy;
import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import Math.EValueCalculator;
import Reports.HistogramVisualizer;
import Utilities.U;

/**
 * Holds peaks and meta-data for a MS/MS peak file.
 * Has many methods for importing peak files.  Has some analysis
 * methods like finding the average peak intensity, etc.
 * @author Brian Risk
 *
 */
public class Spectrum implements Comparable<Spectrum>{

	private ArrayList<Peak> peaks;
	private double maxMass;
	private double precursorMass;
	private double precursorMZ;
	private int id;
	private int charge = 0;
	private File file;
	private String MD5 = null;
	
	//E-Values
	private EValueCalculator eValueCalculator;
	
	//in reality these values should never be less than a positive number
	private double averageIntensity = -1; 
	private double maxIntensity = -1; 
	private double medianIntensity = -1; 
	private double intensity25Percent = -1;
	private double intensity12Percent = -1;
	private double intensity06Percent = -1;
	private double minimumIntensity = -1; 
	private double coverage = -1;
	
	public Spectrum() {
		peaks = new ArrayList<Peak>();
	}
	
	public Spectrum (String s) {
		this(new File(s));
	}
	
	public Spectrum(File file) {
		this(file, false);
	}
	
	public Spectrum(File file, boolean normalizePeaks) {
		this();
		this.file = file;
		Spectrum spectrum = loadSpectra(file).get(0); 
		peaks = spectrum.getPeaks();
		precursorMass = spectrum.getPrecursorMass();
		if (normalizePeaks) normalizePeaks();
	}
	
	public Spectrum(String line, File file) {
		this();
		this.file = file;
		addPrecursorFromString(line);
	}
	
	public void addPeakFromString(String s) {
		try {
			Peak p = new Peak(s);
			if (p.getMass() <= precursorMass) peaks.add(new Peak(s));
		} catch (Exception e) {
			//don't add it if it's bad!
		}	
	}
	
	public void addPrecursorFromString(String s) {
		String [] chunks;
		chunks = s.split("\\s+"); //split on all white space
		precursorMZ = Double.parseDouble(chunks[0]);
		precursorMass = precursorMZ - Definitions.HYDROGEN_MONO;
		if (file.getName().endsWith(".dta")) {
			charge = Integer.parseInt(chunks[1]);
		}
		//assumes txt files are really pkl
		if (file.getName().endsWith(".pkl") || file.getName().endsWith(".txt")) {
			charge = Integer.parseInt(chunks[2]);
			precursorMass *= charge;
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
	
	
	public void sortByIntensity() {
		Peak p;
		for (int i = 0; i < peaks.size(); i++) {
			p = (Peak) peaks.get(i);
			p.setCompareByIntensity();
		}
		Collections.sort(peaks);
	}
	
	public void sortByMass() {
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
		sortByIntensity();
		medianIntensity = getPeak(peaks.size() / 2).getIntensity();
		intensity25Percent = getPeak((int) (peaks.size() * .75)).getIntensity();
		intensity12Percent = getPeak((int) (peaks.size() * .875)).getIntensity();
		intensity06Percent = getPeak((int) (peaks.size() * .9375)).getIntensity();
		sortByMass();
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
			double windowSize = 2 * Properties.peakDifferenceThreshold;
			for (Peak peak: peaks) {
				if (peak.getMass() < precursorMass) {
					upperBound = peak.getMass() + Properties.peakDifferenceThreshold;
					lowerBound = peak.getMass() - Properties.peakDifferenceThreshold;
					if (lowerBound < lastUpperBound) {
						covered += upperBound - lastUpperBound;
					} else {
						covered += windowSize;
					}
					lastUpperBound = upperBound;
				}
			}
			coverage = covered / precursorMass;
		}
		return coverage;
	}

	public Peak getPeak(int i) {return (Peak) peaks.get(i);}
	
	public ArrayList<Peak> getPeaks() {return peaks;}
	
	public double getPrecursorMass() {return precursorMass;}
	public double getPrecursorMZ() {return precursorMZ;}
	
	public int getPeakCount() {return peaks.size();}


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

	public void setAllPeaksToColor(Color c) {
		for (int i = 0; i < getPeakCount(); i++) {getPeak(i).setColor(c);}
	}

	private void cleanPeaks() {
		//before we mess with the peak data, let's make sure we have the MD5
		MD5 = getMD5();
		
		cleanWithWindow();
		
		if (Properties.highIntensityCleaning) 
			cleanPeaksKeepingHighIntensity();
		
		//take exponent of peak intensities
//		if (Properties.defaultScore == Properties.DEFAULT_SCORE_TANDEM_FIT) {
//			for (Peak peak: peaks) {
//				peak.setIntensity(Math.pow(peak.getIntensity(), Properties.peakIntensityExponent));
//			}
//		}
	}
	
	/**
	 * There may be some peaks that don't ever get considered with
	 * scoring mechanisms like tandemFit.  For example, if a very strong peak
	 * is at 1000 Da and another, weaker peak is at 1000.1 Da with nothing close after it
	 * then it will never be used as it will always be overshadowed.  This method
	 * removes those vestigal peaks.
	 */
	private void cleanWithWindow() {
		ArrayList<Peak> retPeaks = new ArrayList<Peak>();
		//add the first peak
		retPeaks.add(peaks.get(0));
		boolean left, right;
		for (int i = 1; i < peaks.size() - 1; i++) {
			left = false;
			right = false;
			for (int j = i - 1; j >=0; j--) {
				if (peaks.get(i).getMass() - peaks.get(j).getMass() <= Properties.peakDifferenceThreshold) {
					if (peaks.get(i).getIntensity() < peaks.get(j).getIntensity()) {
						left = true;
						break;
					}
				} else {
					break;
				}
			}
			if (left) {
				for (int j = i + 1; j < peaks.size(); j++) {
					if ( peaks.get(j).getMass() - peaks.get(i).getMass() <= Properties.peakDifferenceThreshold) {
						if (peaks.get(i).getIntensity() < peaks.get(j).getIntensity()) {
							right = true;
							break;
						}
					} else {
						break;
					}
				}
			}
			if (!(left && right)) {
				retPeaks.add(peaks.get(i));
			}
		}
		//add the final peak
		retPeaks.add(peaks.get(peaks.size() - 1));
		peaks = retPeaks;
	}
	
	//keep the 100 most intense peaks.  Intense!
	private void cleanPeaksKeepingHighIntensity() {
		if (peaks.size() <= Properties.numberOfHighIntensityPeaksToRetain) return;
		sortByIntensity();
		ArrayList<Peak> retPeaks = new ArrayList<Peak>();
		for (int i = peaks.size() - 1; i >= peaks.size() - Properties.numberOfHighIntensityPeaksToRetain; i--) {
			retPeaks.add(peaks.get(i));
		}
		
		peaks = retPeaks;	
		sortByMass();
	}
	
	/**
	 * Finds the peak with maximum intensity.  Goes through
	 * each peak and divides intensity by the max intensity.
	 */
	public void normalizePeaks() {
		double maxIntensity = getMaxIntensity();
		for (Peak peak: peaks) {
			peak.setIntensity(peak.getIntensity() * 100.0 / maxIntensity);
		}
	}
	
	public void print() {
		for (int i = 0; i < peaks.size(); i++) {
			Peak p = (Peak) peaks.get(i);
			System.out.println(p.getMass() + "\t" + p.getIntensity());
		}
		System.out.println("");
	}
	
	public static ArrayList<Spectrum> loadSpectra() {
		ArrayList<Spectrum> spectra = new ArrayList<Spectrum>();
		if (Properties.spectraDirectoryOrFile.isFile()) {
			spectra =  loadSpectra(Properties.spectraDirectoryOrFile);
		} else {
			loadSpectraFromFolder(Properties.spectraDirectoryOrFile, spectra );
		}
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
//						spectrum.sortByMass();
						spectrum.cleanPeaks();
						spectra.add(spectrum);
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
//				spectrum.sortByMass();
				spectrum.cleanPeaks();
				spectra.add(spectrum);
			}
			inBR.close();
		}
		catch (IOException fnfe) {fnfe.printStackTrace(); U.p(inFile.getAbsolutePath());}
		catch (Exception e) {U.p(inFile.getName()); e.printStackTrace(); System.exit(1);}
		return spectra;
	}
	
	/**
	 * ionMathTally is how many ions match with a theoretical peptide.
	 * this uses a stochastic approach to find an upper bound for what the
	 * maximum value of ionMatchTally peak intensities is.
	 * @param ionMatchTally
	 * @return
	 */
	public double getMaxValueForCombination(int ionMatchTally) {
		if (maxValueForCombination[ionMatchTally] == 0) {
//			Random random = new Random();
//			for (int i = 0; i < 10000; i++) {
//				double intensityTotal = 0;
//				for (int j = 0; j < ionMatchTally; j++) {
//					intensityTotal += peaks.get(random.nextInt(peaks.size())).getIntensity();
//				}
//				if (intensityTotal > maxValueForCombination[ionMatchTally]) maxValueForCombination[ionMatchTally] = intensityTotal;
//			}
			sortByIntensity();
			for (int i = peaks.size() - ionMatchTally; i < peaks.size(); i++) {
				if (i < 0) continue;
				maxValueForCombination[ionMatchTally] += peaks.get(i).getIntensity();
			}
			sortByMass();
		}
		return maxValueForCombination[ionMatchTally];
	}
	private double maxValueForCombination[] = new double[100];
	
	public static void main(String args[]) {
		U.p("printing spectrum intensity histograms");
		final int numberOfHistogramBars = 100;
		File dir = new File("SpectraHistograms");
		dir.mkdir();
		ArrayList<Spectrum> spectra = new ArrayList<Spectrum>();
		loadSpectraFromFolder(new File("/Users/risk2/PeppyOverflow/tests/USP/spectra/"), spectra);
		int ionMatchTally = 20;
		Random random = new Random();
		for (int spectrumCount = 0; spectrumCount < 20; spectrumCount++) {
			Spectrum spectrum = spectra.get(spectrumCount);
			spectrum.sortByIntensity();
			double [] histogram = new double[numberOfHistogramBars];
			ArrayList<Peak> peaks = spectrum.getPeaks();
			double min = ionMatchTally * peaks.get(0).getIntensity();
//			double max = ionMatchTally * peaks.get(peaks.size() - 1).getIntensity();
			double max = spectrum.getMaxValueForCombination(ionMatchTally);
			double barWidth = (max - min) / numberOfHistogramBars;
			int bin;
			double intensityTotal;
//			double histoScale = Math.sqrt(peaks.get(peaks.size() - 1).getIntensity() 
//			max = 0;
//			for (int i = 0; i < 200000; i++) {
//				intensityTotal = 0;
//				for (int j = 0; j < ionMatchTally; j++) {
//					intensityTotal += peaks.get(random.nextInt(peaks.size())).getIntensity();
//				}
//				if (intensityTotal > max) max = intensityTotal;
//			}
//			barWidth = (max - min) / numberOfHistogramBars;
			
			for (int i = 0; i < 100000; i++) {
				intensityTotal = 0;
				for (int j = 0; j < ionMatchTally; j++) {
					intensityTotal += peaks.get(random.nextInt(peaks.size())).getIntensity();
				}
				bin = (int) Math.floor(( intensityTotal - min) / barWidth);
				if (bin < numberOfHistogramBars) {
					histogram[bin]++;
				} else {
					histogram[numberOfHistogramBars - 1]++;
				}
			}
//			for (Peak peak: peaks) {
//				bin = (int) Math.floor((transform(peak.getIntensity()) - min) / barWidth);
//				if (bin < numberOfHistogramBars) {
//					histogram[bin]++;
//				} else {
//					histogram[numberOfHistogramBars - 1]++;
//				}
//			}
//			//probabilities
//			for (int i=0; i<numberOfHistogramBars; i++) {
//				histogram[i] /= peaks.size();
//			}
			//make survival
			for (int i=numberOfHistogramBars - 2; i >=0 ; i--) {
				histogram[i] += histogram[i + 1];
			}
			//log of the bars
			for (int i=0; i<numberOfHistogramBars; i++) {
				if (histogram[i] > 0) {
					histogram[i] = Math.abs(Math.log(histogram[i]));
				}
			}
			try {
				String fileName = U.getFileNameWithoutSuffix(spectrum.getFile()) + ".jpg";
				HistogramVisualizer.drawHistogram(histogram, 600, 800, new File(dir, fileName));
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		U.p("done");
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
	

	public int compareTo(Spectrum spectrum) {
		if (precursorMass < spectrum.getPrecursorMass()) return -1;
		if (precursorMass > spectrum.getPrecursorMass())  return 1;
		return  0;
	}
	
	public void clearEValues() {
		eValueCalculator = null;
	}
	
	public void calculateEValues(ArrayList<Match> matches, ArrayList<Match> topMatches) {
		if (eValueCalculator == null) {
			eValueCalculator = new EValueCalculator(matches, topMatches);
		} else {
			eValueCalculator.addScores(matches, topMatches);
		}
	}
	
	public double getEValue(double score) {
		return eValueCalculator.getEValue(score);
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
	
	

}

