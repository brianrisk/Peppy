package Peppy;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Scanner;

import Math.MassError;

import uk.ac.ebi.jmzml.model.mzml.BinaryDataArray;
import uk.ac.ebi.jmzml.model.mzml.BinaryDataArrayList;
import uk.ac.ebi.jmzml.model.mzml.CVParam;
import uk.ac.ebi.jmzml.model.mzml.ParamGroup;
import uk.ac.ebi.jmzml.model.mzml.Precursor;
import uk.ac.ebi.jmzml.model.mzml.PrecursorList;
import uk.ac.ebi.jmzml.model.mzml.SelectedIonList;
import uk.ac.ebi.jmzml.xml.io.MzMLObjectIterator;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;

public class SpectrumLoader {

	/*Public Static methods for loading spectra*/
	/**
	 * loadSpectra loads in spectra from the established value in the Properties file.  It returns an ArrayList of Spectrum that is the same size (minus removed peaks)
	 * as the number of spectrum in the properties file.  It removes any spectra below the properties file established value of minimumNumberOfPeaksForAValidSpectrum.
	 * 
	 * This methods defaults to the spectra File set by properties
	 * 
	 * @return 
	 */
	public static ArrayList<Spectrum> loadSpectra() {

		ArrayList<Spectrum> spectra = new ArrayList<Spectrum>();
		
		//Load spectra based on whether it is a single file or a directory
		if (Properties.spectraDirectoryOrFile.isFile()) {
			spectra =  loadSpectra(Properties.spectraDirectoryOrFile);
		} else {
			loadSpectraFromFolder(Properties.spectraDirectoryOrFile, spectra );
		}//else
		
		/* remove spectra with small amount of peaks */
		int removeTally = 0;
		for (int i = 0; i < spectra.size(); i++) {
			if (spectra.get(i).getPeakCount() < Properties.minimumNumberOfPeaksForAValidSpectrum) {
				//Remove a peak that has to few peaks
				spectra.remove(i);
				//Compensate for the newly removed spectra
				i--;
				//Keep track of how many are moved
				removeTally++;
			}//if
		}//else
		
		//Display how many peaks are removed
		if (removeTally > 0) {
			U.p("Removed " + removeTally + " spectra with less than " + Properties.minimumNumberOfPeaksForAValidSpectrum + " peaks.");
		}//if
		
		/* set spectra id */
		for (int i = 0; i < spectra.size(); i++) {
			spectra.get(i).setId(i);
		}//for
		
		
		return spectra;
	
	}//loadSpectra
	
	/**
	 * A utility method to load in a DTA/PK file.  This is a convenience method for ease of programming. It calls loadSpectra(File inFile).
	 * @param fileName  A string representing the locaiton of the .dta file to load
	 * @return An ArrayList of all the spectra
	 */
	public static ArrayList<Spectrum> loadSpectra(String fileName) {
		File inFile = new File(fileName);
		return loadSpectra(inFile);
	}//loadSpectra
	
	/**
	 * A utility method to load in a DTA/PKL file.  There may be more than
	 * one spectrum in a file, so this method returns an ArrayList.
	 * 
	 *This is the workhorse method for loading spectra.  In the end every spectra is loaded using this method
	 * @param inFile
	 * @return A ArrayList of all the spectra
	 */
	public static ArrayList<Spectrum> loadSpectra(File inFile) {
		//if the spectra are in a .mzml file, then use a different method and return;
		if(inFile.getAbsolutePath().toLowerCase().endsWith(".mzml") ){
			return loadMZMLSpectra(inFile);
		}//if
		
		if(inFile.getAbsolutePath().toLowerCase().endsWith(".mgf")){
			return loadMGFSpectra(inFile);
		}
		
		if(inFile.getAbsolutePath().toLowerCase().endsWith(".dta")){
			return loadDTASpectra(inFile);
		}
		
		if(inFile.getAbsolutePath().toLowerCase().endsWith(".pkl") || inFile.getAbsolutePath().toLowerCase().endsWith(".txt")){
			return loadPKLSpectra(inFile);
		}
		
		
		
		return null;
	}//loadSpectra
	
	private static ArrayList<Spectrum> loadPKLSpectra(File inFile){
		ArrayList<Spectrum> out = new ArrayList<Spectrum>();

		ArrayList<Peak> peaks = new ArrayList<Peak>();
			try {
				//BufferedReader to upload the spectra file
				BufferedReader inBR = new BufferedReader(new FileReader(inFile));
				//Temporary variable to store a single line
				String line = inBR.readLine();
				//Value to store 
				boolean fileOpen = false;
				//Single spectrum object to store information before putting it into the spectra output variable
				Spectrum spectrum = null;
				
				//Iterate through the file and grab each spectrum
				while (line != null) {
					//The end of a spectrum/file is found
					if (line.trim().equals("")) {
						
						//This is the end of the file, so close it out and get the last spectrum cleaned up
						if (fileOpen) {
							fileOpen = false;
							/* clean peaks also sorts the peaks by mass */
							
							spectrum.setPeaks(peaks);
							cleanPeaks(spectrum);
							spectrum.setTitle(inFile.getAbsolutePath().substring(inFile.getAbsolutePath().lastIndexOf('/')));
							out.add(spectrum);
						}
						//An end of a line has not been reached
					}else{
						
						//If the file is already open, start adding lines from the file
						if (fileOpen) {
							try {
								Peak p = new Peak(line);
								peaks.add(p);
							} catch (Exception e) {
								//don't add it if it's bad!
							}
							
							//else this is the first line of a new spectrum
						} else {
							fileOpen = true;
							spectrum = new Spectrum();
							spectrum.setFile(inFile);
							
							//Add in precursor information
							String [] chunks = line.split("\\s+"); //split on all white space
							spectrum.setPrecursorMZ(Double.parseDouble(chunks[0]));
							spectrum.setMass(spectrum.getPrecursorMZ() - Definitions.HYDROGEN_MONO);
							spectrum.setCharge(Integer.parseInt(chunks[2]));
							spectrum.setMass(spectrum.getMass() * spectrum.getCharge());
							
							//Reset the peaks for this spectra
							peaks = new ArrayList<Peak>();

						}//else

					}//Middle else
					line = inBR.readLine();
				}//while loop
				
				//If the file is still open at this point it must be the last spectrum since the lines in the file have run out
				if (fileOpen) {
					/* clean peaks also sorts the peaks by mass */
					if (spectrum.isValid()) {
						spectrum.setPeaks(peaks);
						cleanPeaks(spectrum);
						spectrum.setTitle(inFile.getAbsolutePath().substring(inFile.getAbsolutePath().lastIndexOf('/')));
						out.add(spectrum);
					}//if
					
				}//if file is open
				
				//Close out the buffered reader
				inBR.close();
				
			}//try
			catch (IOException fnfe) {fnfe.printStackTrace(); U.p(inFile.getAbsolutePath());}
			catch (Exception e) {U.p(inFile.getName()); e.printStackTrace(); System.exit(1);}
		
		return out;
	}//loadPKLSpectra
	
	private static ArrayList<Spectrum> loadDTASpectra(File inFile){
		ArrayList<Spectrum> out = new ArrayList<Spectrum>();

		ArrayList<Peak> peaks = new ArrayList<Peak>();
			try {
				//BufferedReader to upload the spectra file
				BufferedReader inBR = new BufferedReader(new FileReader(inFile));
				//Temporary variable to store a single line
				String line = inBR.readLine();
				//Value to store 
				boolean fileOpen = false;
				//Single spectrum object to store information before putting it into the spectra output variable
				Spectrum spectrum = null;
				
				//Iterate through the file and grab each spectrum
				while (line != null) {
					//The end of a spectrum/file is found
					if (line.trim().equals("")) {
						
						//This is the end of the file, so close it out and get the last spectrum cleaned up
						if (fileOpen) {
							fileOpen = false;
							/* clean peaks also sorts the peaks by mass */
							spectrum.setPeaks(peaks);
							cleanPeaks(spectrum);
							spectrum.setTitle(inFile.getAbsolutePath().substring(inFile.getAbsolutePath().lastIndexOf('/')));
							out.add(spectrum);
						}
						//An end of a line has not been reached
					}else{
						
						//If the file is already open, start adding lines from the file
						if (fileOpen) {
							try {
								Peak p = new Peak(line);
								peaks.add(p);
							} catch (Exception e) {
								//don't add it if it's bad!
							}
						} 
						//else this is the first line of a new spectrum
						else {
							fileOpen = true;
							spectrum = new Spectrum();
							spectrum.setFile(inFile);

							//Add in precusor infomration
							String [] chunks = line.split("\\s+"); //split on all white space
							spectrum.setPrecursorMZ(Double.parseDouble(chunks[0]));
							spectrum.setMass(spectrum.getPrecursorMZ() - Definitions.HYDROGEN_MONO);
							spectrum.setCharge(Integer.parseInt(chunks[1]));

							//Reset the peaks for this spectra
							peaks = new ArrayList<Peak>();

						}//else

					}//Middle else
					line = inBR.readLine();
				}//while loop
				
				//If the file is still open at this point it must be the last spectrum since the lines in the file have run out
				if (fileOpen) {
					/* clean peaks also sorts the peaks by mass */
					if (spectrum.isValid()) {
						spectrum.setPeaks(peaks);
						cleanPeaks(spectrum);
						spectrum.setTitle(inFile.getAbsolutePath().substring(inFile.getAbsolutePath().lastIndexOf('/')));
						out.add(spectrum);
					}//if
					
				}//if file is open
				
				//Close out the buffered reader
				inBR.close();
				
			}//try
			catch (IOException fnfe) {fnfe.printStackTrace(); U.p(inFile.getAbsolutePath());}
			catch (Exception e) {U.p(inFile.getName()); e.printStackTrace(); System.exit(1);}


		return out;
	}//loadDTASpectra
	
	private static ArrayList<Spectrum> loadMZMLSpectra(File inFile){
//		U.p("Loading up mzml spectra!");
		ArrayList<Spectrum> out = new ArrayList<Spectrum>();
		
		MzMLUnmarshaller unmarshaller = new MzMLUnmarshaller(inFile);
		
		
		MzMLObjectIterator<uk.ac.ebi.jmzml.model.mzml.Spectrum> spectrumIterator = unmarshaller.unmarshalCollectionFromXpath("/run/spectrumList/spectrum", uk.ac.ebi.jmzml.model.mzml.Spectrum.class);
		while (spectrumIterator.hasNext()){
		  //read next spectrum from XML file
			uk.ac.ebi.jmzml.model.mzml.Spectrum spectrum = spectrumIterator.next();

		//Create a list to retreive mass and intensity data
		  BinaryDataArrayList bdal = spectrum.getBinaryDataArrayList();
		  List<BinaryDataArray> lbd = bdal.getBinaryDataArray();
		  
		  //If a spectra has both mass and intensity list add it to the output
		  if(lbd.size() >= 2){
			  
			  //Create mass and Inetesity list to access for the creation of this spectrum
			  BinaryDataArray mass = lbd.get(0);
			  BinaryDataArray intensity = lbd.get(1);
			  //Create a dummy spectrum
			  Spectrum temp = new Spectrum();
			  temp.setTitle(inFile.getAbsolutePath().substring(inFile.getAbsolutePath().lastIndexOf('/')));
			  temp.setFile(inFile);
			  
			  //Add in precursor
				PrecursorList pcl = spectrum.getPrecursorList();
				
				//Get the precursor information if this spectrum has it
				//Should this be modified to require precursor information in a spectrum.  right now a spectra cannot have it and wills till be created.
				if(pcl != null){
					List<Precursor> lpc = pcl.getPrecursor();
					SelectedIonList sil = lpc.get(0).getSelectedIonList();
					List<ParamGroup> lsi = sil.getSelectedIon();
					List<CVParam> lcv = lsi.get(0).getCvParam();
					//This is the same as .PKL format files
					//0: This grabs selected ion m/z
					//1: This grabs peak intensity
					//2: This grabs charge state
					temp.setPrecursorMZ(Double.parseDouble(lcv.get(0).getValue()));
					temp.setMass(temp.getPrecursorMZ() - Definitions.HYDROGEN_MONO);
					temp.setCharge(Integer.parseInt(lcv.get(2).getValue()));
					
				}//if pcl is not null(The spectra has precursor information
			
				ArrayList<Peak> peaks = new ArrayList<Peak>();
			  //Add in mass and Intensity
			  for(int j = 0; j < mass.getBinaryDataAsNumberArray().length;j++){
				  
				  //Store the mass + intensity of this spectrum
	
					 
				 try {
					Peak p = new Peak(mass.getBinaryDataAsNumberArray()[j] + "\t" + intensity.getBinaryDataAsNumberArray()[j]);
					peaks.add(p);
				} catch (Exception e) {
					//don't add it if it's bad!
				}



			  }//for
			  
			  //Save the peaks read int from the file
			  temp.setPeaks(peaks);
			  
			  //Store the spectrum into output
			  out.add(temp);
			  
		  }//if size == 2
		  
		}//While iterator hasNext
		
		return out;
	}//loadMZML Spectra
	
	private static ArrayList<Spectrum> loadMGFSpectra(File inFile){
		ArrayList<Spectrum> out = new ArrayList<Spectrum>();
//		U.p("Uploading a mgf file");
		Spectrum temp = new Spectrum();
		
		
		//Prase the file and create a series of spectrumt o add to the output
		try {
			Scanner s = new Scanner(inFile);
			
			while(s.hasNextLine()){
				String line = s.nextLine();
				
				//Skip comment lines
				if(line.startsWith("#")){
					continue;
				}//if comment line
				//SKip empty lines with just newline characters
				if(line.length() <= 1){
					continue;
				}
				
				if(line.equals("BEGIN IONS")){
					line = s.nextLine();
					//Loop through and load up a spectra
					
					String[] titleLine = line.split("=");
					String[] brokenTitleLine = titleLine[1].split(" ");
					
					
					/*Parse the header lines from this spectra*/
					String title = brokenTitleLine[0].substring(brokenTitleLine[0].indexOf(':')  + 1, brokenTitleLine[0].length());

					String scans = s.nextLine().split("=")[1];
					String RTINSECONDS = s.nextLine().split("=")[1];
					
					
					String charge = s.nextLine().split("=")[1];
					String pepmass = s.nextLine().split("=")[1];
					
					temp.setTitle(title);
					temp.setFile(inFile);
					
					//Create a string based on the mass of this peptide

					charge = charge.substring(0, charge.length() - 1);
					//Add in precursor
					
					temp.setPrecursorMZ(Double.parseDouble(pepmass));
					temp.setMass(temp.getPrecursorMZ() - Definitions.HYDROGEN_MONO);
					temp.setCharge(Integer.parseInt(charge));
					

					
					ArrayList<Peak> peaks = new ArrayList<Peak>();
					while(s.hasNextLine()){
						line = s.nextLine();
						if(line.equals("END IONS")){
							break;
						}
						String [] chunks = line.split(" |\t");
						
						//Only add value peaks that have a positive mass
				
						try {
							Peak p = new Peak(line);
							peaks.add(p);
						} catch (Exception e) {
							//don't add it if it's bad!
						}	

						
					}//while
					
					//Set the peaks to this spectrum object
					temp.setPeaks(peaks);
					
					//Add this spectrum to the output
					out.add(temp);
					
					temp = new Spectrum();
					
					
				}//BEGIN IONS
				
			}//while
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		
		return out;
	}//loadMGFSpectaInFile
	

	/**
	 * loadSpectraFromFolder takes the following steps:
	 * 1) Recursively goes through folder.
	 * 2) Finds all files that end in .dta or .pkl or .txt
	 * 3) Extracts the spectra from those files
	 * 4) Adds all those spectra to one big ArrayList and returns that.
	 * @param folder  A file object that represents the folder containing spectra for this method to use
	 * @param spectra  The ArrayList to add spectra to
	 */
	public static void loadSpectraFromFolder(File folder, ArrayList<Spectrum> spectra) { 
		File [] files = folder.listFiles();
		
		//Loop through and grab all spectra from the folder
		for (int i = 0; i < files.length; i++) {
			if (files[i].isHidden()) continue;
			if (files[i].isDirectory()) {
				loadSpectraFromFolder(files[i], spectra);
				continue;
			}
			String fileName = files[i].getName().toLowerCase();
			if (fileName.toLowerCase().endsWith(".dta") || fileName.toLowerCase().endsWith(".pkl") || fileName.toLowerCase().endsWith(".txt") || 
					fileName.toLowerCase().endsWith(".mzml") || fileName.toLowerCase().endsWith(".mgf")) {
				spectra.addAll(loadSpectra(files[i]));
			}
		}
		//giving all spectra an ID
		for (int i = 0; i < spectra.size(); i++) {
			spectra.get(i).setId(i);
		}
	}//loadSpectraFromFOlder
	
	/**
	 * This is slightly more general than the method it overloads.  It creates 
	 * a spectrum ArrayList which it passes to the other loadSpectraFromFolder
	 * This method uses a string instead of a File object as input
	 * @param fileName
	 * @return 
	 */
	public static ArrayList<Spectrum> loadSpectraFromFolder(String fileName) {
		File inFile = new File(fileName);
		return loadSpectraFromFolder(inFile);
	}//loadSpectraFromFolder
	
	
	/**
	 * This is slightly more general than the method it overloads.  It creates 
	 * a spectrum ArrayList which it passes to the other loadSpectraFromFolder
	 * @param fileName
	 * @return
	 */
	public static ArrayList<Spectrum> loadSpectraFromFolder(File inFile) {
		ArrayList<Spectrum> spectra = new ArrayList<Spectrum>();
		loadSpectraFromFolder(inFile, spectra);
		return spectra;
	}//loadSpectraFromFolder
	

	/**
	 * '
	 * This method is .
	 * recursively goes through folder.
	 * finds all files that end in .dta or .pkl or .mzml or .mgf
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
			if (fileName.toLowerCase().endsWith(".dta") || fileName.toLowerCase().endsWith(".pkl") || fileName.toLowerCase().endsWith(".txt") || 
					fileName.toLowerCase().endsWith(".mzml") || fileName.toLowerCase().endsWith(".mgf")) {
				spectraFiles.add(files[i]);
			}
		}
	}
	
	
/*End Public Static methods for loading spectra*/
	
	
	/**
	 * cleanPeaks ensures a spectrum has a MD5, and then keeps the strongest peak in each region and sort the peaks by mass.
	 */
	public static void cleanPeaks(Spectrum s) {
		
		//before we mess with the peak data, let's make sure we have the MD5
		s.setMD5(s.getMD5());
	

		s.setPeaks(keepStrongestPeakInRegions(s));
		sortPeaksByMass(s);
		
	}//cleanPeaks
	

	
	/**
	 * keepStrongestPeakInRegions sorts the peaks by intensity, and then ensures that know two peaks are within each other by a factor of the fragmentTolerence
	 * variable from the Properties file.  
	 */
	private static ArrayList<Peak> keepStrongestPeakInRegions(Spectrum s) {
		ArrayList<Peak> peaks = s.getPeaks();
		sortPeaksByIntensity(s);
		int start = s.getPeaks().size() - 1;
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
					}//mass < upperBound
				}//mass > lower bound
			}//for start to stop
		}//for start ot stop
		
		return peaks;
	}//keepStrongestPeakInRegions
	
	
/*Sorting methods*/
	
	public static void sortPeaksByIntensity(Spectrum s) {
		ArrayList<Peak> peaks = s.getPeaks();
		Peak p;
		for (int i = 0; i < peaks.size(); i++) {
			p = (Peak) peaks.get(i);
			p.setCompareByIntensity();
		}
		Collections.sort(peaks);
		s.setPeaks(peaks);
	}//sortPeaksByIntensity
	
	public static void sortPeaksByMass(Spectrum s) {
		ArrayList<Peak> peaks = s.getPeaks();
		Peak p;
		for (int i = 0; i < peaks.size(); i++) {
			p = (Peak) peaks.get(i);
			p.setCompareByMass();
		}
		Collections.sort(s.getPeaks());
		s.setPeaks(peaks);
	}//sortPeaksByMass

/*End Sorting methods*/
	
}
