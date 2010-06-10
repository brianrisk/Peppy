package Peppy;
import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collections;
import java.util.ArrayList;

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
	private int id;
	private int charge = 0;
	private File file;
	
	//in reality these values should never be less than a positive number
	private double averageIntensity = -1; 
	private double maxIntensity = -1; 
	private double medianIntensity = -1; 
	private double minimumIntensity = -1; 
	
	public Spectrum() {
		peaks = new ArrayList<Peak>();
	}
	
	public Spectrum(File file) {
		this();
		this.file = file;
		Spectrum spectrum = loadSpectra(file).get(0); 
		peaks = spectrum.getPeaks();
		precursorMass = spectrum.getPrecursorMass();
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
		} catch (MalformedPeakException mpe) {
			//don't add it if it's bad!
		}	
	}
	
	public void addPrecursorFromString(String s) {
		String [] chunks;
		chunks = s.split("\\s+"); //split on all white space
		precursorMass = Double.parseDouble(chunks[0]);
		if (file.getName().endsWith(".dta")) {
			charge = Integer.parseInt(chunks[1]);
			precursorMass -= Definitions.HYDROGEN_MONO;
		}
		//assumes txt files are really pkl
		if (file.getName().endsWith(".pkl") || file.getName().endsWith(".txt")) {
			charge = Integer.parseInt(chunks[2]);
			precursorMass = charge * (precursorMass - Definitions.HYDROGEN_MONO);
		}
	}
	
	/* returns a double, but we're really interested in setting the max */
	public double getMaxMass() {
		sortByMass();
		Peak p = (Peak) peaks.get(peaks.size() - 1);
		maxMass = p.getMass();
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
	
	//returns the mass for peak i
	public double getMass(int i) {
		Peak p = (Peak) peaks.get(i);
		return p.getMass();
	}
	
	public double calculateAverageIntensity() {
		double out = 0.0;
		Peak p;
		for (int i = 0; i < peaks.size(); i++) {
			p = (Peak) peaks.get(i);
			out += p.getIntensity();
		}
		out /= peaks.size();
		averageIntensity = out;
		return averageIntensity;
	}
	
	public double calculateMaxIntensity() {
		sortByIntensity();
		Peak p = getPeak(getPeakCount() - 1);
		sortByMass();
		maxIntensity = p.getIntensity();
		return maxIntensity;
	}
		
	public double calculateMedianIntensity() {
		sortByIntensity();
		Peak p = getPeak(peaks.size() / 2);
		sortByMass();
		medianIntensity = p.getIntensity();
		return medianIntensity;
	}
	
	public double calculateMinimumIntensity() {
		sortByIntensity();
		Peak p = getPeak(0);
		sortByMass();
		minimumIntensity = p.getIntensity();
		return minimumIntensity;
	}
	
	/**
	 * If the average has not been calculated already, it calculates it
	 * @return
	 */
	public double getCalculatedAverageIntensity() {
		if (averageIntensity < 0) return calculateAverageIntensity();
		return averageIntensity;
	}
	
	/**
	 * If the maximum has not been found already, it finds it
	 * @return
	 */
	public double getCalculatedMaxIntensity() {
		if (maxIntensity < 0) return calculateMaxIntensity();
		return maxIntensity;
	}
	
	/**
	 * If the minimum has not been found already, it finds it
	 * @return
	 */
	public double getCalculatedMinimumIntensity() {
		if (minimumIntensity < 0) return calculateMinimumIntensity();
		return minimumIntensity;
	}
	
	/**
	 * If the median has not been found already, it finds it
	 * @return
	 */
	public double getCalculatedMedianIntensity() {
		if (medianIntensity < 0) return calculateMedianIntensity();
		return medianIntensity;
	}
	
	
	public Peak getPeak(int i) {return (Peak) peaks.get(i);}
	
	public ArrayList<Peak> getPeaks() {return peaks;}
	
	public double getPrecursorMass() {return precursorMass;}
	
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

	public void setAllPeaksToColor(Color c) {
		for (int i = 0; i < getPeakCount(); i++) {getPeak(i).setColor(c);}
	}

	public void cleanPeaks() {
		if (Properties.highIntensityCleaning) 
			cleanPeaksKeepingHighIntensity();
		if (Properties.localMaximaCleaning) 
			cleanPeaksKeepingLocalMaxima();
	}
	
	private void cleanPeaksKeepingLocalMaxima() {
		ArrayList<Peak> retPeaks = new ArrayList<Peak>();
		int i, peakCount = peaks.size();
		Peak thisPeak, nextPeak, prevPeak;
		float gap = 5.0f;
		retPeaks.add(peaks.get(0)); //add first peak
		for (i=1; i<peakCount - 1; i++) {
			thisPeak = (Peak) peaks.get(i);
			nextPeak = (Peak) peaks.get(i+1);
			prevPeak = (Peak) peaks.get(i-1);

			//add the peak if it is on the edge of a gap
			if	(
				(thisPeak.getIntensity() - prevPeak.getIntensity() > gap) ||
				(nextPeak.getIntensity()  - thisPeak.getIntensity()  > gap)
			){
				retPeaks.add(thisPeak);
				continue;
			}
			if ( (thisPeak.getIntensity() > nextPeak.getIntensity()) && (thisPeak.getIntensity() > prevPeak.getIntensity()) ) {
				retPeaks.add(thisPeak);
			}
		}
		if (peakCount > 1) {
			retPeaks.add(peaks.get(peakCount - 1));
		}
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
	 * a utility method to load in a DTA/PKL file
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
	}
	
	/**
	 * This is slightly more general than the method it overloads.  It creates 
	 * a spectrum ArrayList which it passes to the other loadSpectraFromFolder
	 * @param fileName
	 * @return
	 */
	public static ArrayList<Spectrum> loadSpectraFromFolder(String fileName) {
		ArrayList<Spectrum> spectra = new ArrayList<Spectrum>();
		File inFile = new File(fileName);
		loadSpectraFromFolder(inFile, spectra);
		return spectra;
	}
	

	public int compareTo(Spectrum spectrum) {
		if (precursorMass < spectrum.getPrecursorMass()) return -1;
		if (precursorMass > spectrum.getPrecursorMass())  return 1;
		return  0;
	}



}

