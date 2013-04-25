package Peppy;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Scanner;

import org.apache.log4j.Level;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;

import uk.ac.ebi.jmzml.model.mzml.BinaryDataArray;
import uk.ac.ebi.jmzml.model.mzml.BinaryDataArrayList;
import uk.ac.ebi.jmzml.model.mzml.CVParam;
import uk.ac.ebi.jmzml.model.mzml.MzML;
import uk.ac.ebi.jmzml.model.mzml.ParamGroup;
import uk.ac.ebi.jmzml.model.mzml.Precursor;
import uk.ac.ebi.jmzml.model.mzml.PrecursorList;
import uk.ac.ebi.jmzml.model.mzml.SelectedIonList;
import uk.ac.ebi.jmzml.xml.io.MzMLObjectIterator;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.JMzReaderException;
import uk.ac.ebi.pride.tools.ms2_parser.Ms2File;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLParsingException;
import Math.MassError;


/**
 * 
 * @author Brian Risk 
 *
 */
public class SpectrumLoader {

	
	//A list of accepted file extensions peppy takes in.  When this is updated don't forget to update loadSpectra
	private static String[] fileExtensions = {".dta", ".pkl", ".mzml.xml", ".mzml", ".mgf", ".mz.xml", ".mzxml", ".mzxml.xml", ".ms2"};
	
	
	/*Public Static methods for loading spectra*/
	/**
	 * loadSpectra loads in spectra from the established value in the Properties file.  It returns an ArrayList of Spectrum that is the same size (minus removed peaks)
	 * as the number of spectrum in the properties file.  It removes any spectra below the properties file established value of minimumNumberOfPeaksForAValidSpectrum.
	 * 
	 * This methods defaults to the spectra File set by properties
	 * 
	 * @return An ArrayList of Spectrum objects to loadSpectra
	 */
	public static ArrayList<Spectrum> loadSpectra() {
		return loadSpectra(Properties.spectraDirectoryOrFile);
	
	}//loadSpectra
	
	/**
	 * A utility method to load in a DTA/PK file.  This is a convenience method for ease of programming. It calls loadSpectra(File inFile).
	 * @param fileName  A string representing the location of the .dta file to load
	 * @return An ArrayList of all the spectra
	 */
	public static ArrayList<Spectrum> loadSpectra(String fileName) {
		File inFile = new File(fileName);
		return loadSpectra(inFile);
	}//loadSpectra
	
	
	/**
	 * There may be more than
	 * one spectrum in a file, so this method returns an ArrayList.
	 * 
	 *This is the workhorse method for loading spectra.  In the end every spectrum is loaded using this method
	 * @param inFile
	 * @return An ArrayList of all the spectra
	 */
	public static ArrayList<Spectrum> loadSpectra(File inFile) {
		ArrayList<Spectrum> spectra = new ArrayList<Spectrum>();
		
		if (inFile.isFile()) {
			spectra = loadSpectraFromFile(inFile);
		} else {
			spectra = loadSpectraFromFolder(inFile);
		}
		
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
		
		/*
		 * Giving each spectrum its own ID.
		 * For files with only one spectrum, such as a folder
		 * full of DTA files, this will be a lot of "0" IDs.
		 */
		if (spectra != null) {
			for (int spectrumID = 0; spectrumID < spectra.size(); spectrumID++) {
				spectra.get(spectrumID).setId(spectrumID);
			}
		}
		
		return spectra;

	}//loadSpectra
	
	
	public static ArrayList<Spectrum> loadSpectraFromFile(File inFile) {
		
			ArrayList<Spectrum> spectra = new ArrayList<Spectrum>();
			
		
			//if the spectra are in a .mzml file, then use a different method and return;
			if(inFile.getAbsolutePath().toLowerCase().endsWith(".mzml") || inFile.getAbsolutePath().toLowerCase().endsWith(".mzml.xml") ){
				spectra = loadMZMLSpectra(inFile);
			}//if
			
			if(inFile.getAbsolutePath().toLowerCase().endsWith(".mgf")){
				spectra = loadMGFSpectra(inFile);
			}
			
			if(inFile.getAbsolutePath().toLowerCase().endsWith(".dta")){
				spectra = loadDTASpectra(inFile);
			}

			
			if(inFile.getAbsolutePath().toLowerCase().endsWith(".pkl")){
				spectra = loadPKLSpectra(inFile);
			}
			if(inFile.getAbsolutePath().toLowerCase().endsWith(".mz.xml") || inFile.getAbsolutePath().toLowerCase().endsWith(".mzxml")
						|| inFile.getAbsolutePath().toLowerCase().endsWith(".mzxml.xml")){
				spectra = loadMzXMLSpectra(inFile);
			}
			
			if(inFile.getAbsolutePath().toLowerCase().endsWith(".ms2")){
				spectra = loadMS2Spectra(inFile);
			}
			
			/* for files that have more than one spectrum */
			for (int i = 0; i < spectra.size(); i++) {
				spectra.get(i).setId(i);
			}
			
			return spectra;
	}
	
	
	

	
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
						spectrum.setFile(inFile);
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
							spectrum.setFile(inFile);
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
	
	/**
	 * 
	 * Uses jmzml http://code.google.com/p/jmzml/
	 * @param inFile
	 * @return
	 */
	private static ArrayList<Spectrum> loadMZMLSpectra(File inFile){
		ArrayList<Spectrum> out = new ArrayList<Spectrum>();

		/* turning of logging because it is used extensively in the jmzML code  */
		List<Logger> loggers = Collections.<Logger>list(LogManager.getCurrentLoggers());
		loggers.add(LogManager.getRootLogger());
		for ( Logger logger : loggers ) {
			logger.setLevel(Level.ERROR);
		}

		MzMLUnmarshaller unmarshaller = new MzMLUnmarshaller(inFile);	
//		MzML completeMzML = unmarshaller.unmarshall();
//		List<uk.ac.ebi.jmzml.model.mzml.Spectrum> allSpectra = completeMzML.getRun().getSpectrumList().getSpectrum();
		MzMLObjectIterator<uk.ac.ebi.jmzml.model.mzml.Spectrum> spectrumIterator = unmarshaller.unmarshalCollectionFromXpath("/run/spectrumList/spectrum", uk.ac.ebi.jmzml.model.mzml.Spectrum.class);

		
		System.out.println("Number of spectrum elements: " + unmarshaller.getObjectCountForXpath("/run/spectrumList/spectrum"));

		int spectrumCount = 0;
		long millisecondsThen = System.currentTimeMillis();
		long millisecondsNow;
		int skippedTotal = 0;
		
		uk.ac.ebi.jmzml.model.mzml.Spectrum mzMLSpectrum;
		while (spectrumIterator.hasNext()){
//		for (uk.ac.ebi.jmzml.model.mzml.Spectrum mzMLSpectrum: allSpectra) {

			spectrumCount++;
			if (spectrumCount == 10000) {
//				U.p("skippedTotal: " + skippedTotal);
//				return out;
			}
			if (spectrumCount % 100 == 0) {
				millisecondsNow = System.currentTimeMillis();
//				U.p("spectrumCount: " + spectrumCount + ", " + (millisecondsNow - millisecondsThen));
				millisecondsThen = millisecondsNow;
			}

			/* read next spectrum from XML file */
			mzMLSpectrum = spectrumIterator.next();

			/* Create a list to retrieve mass and intensity data */
			BinaryDataArrayList bdal = mzMLSpectrum.getBinaryDataArrayList();
			
			List<BinaryDataArray> lbd = bdal.getBinaryDataArray();

			/* If a spectra has both mass and intensity list add it to the output */
			if(lbd.size() >= 2){
				
				//Create an empty spectrum
				Spectrum spectrumUnderConstruction = new Spectrum();
				spectrumUnderConstruction.setTitle(inFile.getAbsolutePath().substring(inFile.getAbsolutePath().lastIndexOf('/')));
				spectrumUnderConstruction.setFile(inFile);
				
				/* retrieve precursor values */
				PrecursorList precursorList = mzMLSpectrum.getPrecursorList();
				if(precursorList != null){
					List<Precursor> lpc = precursorList.getPrecursor();
					SelectedIonList sil = lpc.get(0).getSelectedIonList();
					List<ParamGroup> lsi = sil.getSelectedIon();
					List<CVParam> lcv = lsi.get(0).getCvParam();

					for(int i = 0; i < lcv.size(); i++){
						if(lcv.get(i).getName().equals("selected ion m/z")){
							spectrumUnderConstruction.setPrecursorMZ(Double.parseDouble(lcv.get(i).getValue()));
						}
						if(lcv.get(i).getName().equals("charge state")){
							spectrumUnderConstruction.setCharge(Integer.parseInt(lcv.get(i).getValue()));
						}
					}

					/* determining mass and setting that property */
					spectrumUnderConstruction.setMass(spectrumUnderConstruction.getPrecursorMZ() - Definitions.HYDROGEN_MONO);
					spectrumUnderConstruction.setMass(spectrumUnderConstruction.getMass() * spectrumUnderConstruction.getCharge());

				} else {
					skippedTotal++;
					continue;
				}
		
//				/* Create mass and intensity list to access for the creation of this spectrum */
//				BinaryDataArray mass = lbd.get(0);
//				BinaryDataArray intensity = lbd.get(1);
//				Number [] massArray =  mass.getBinaryDataAsNumberArray();
//				Number [] intensityArray =  intensity.getBinaryDataAsNumberArray();
//				ArrayList<Peak> peaks = new ArrayList<Peak>(massArray.length);
//				for(int j = 0; j < massArray.length;j++){
//					Peak p = new Peak(massArray[j].floatValue(), intensityArray[j].floatValue());
//					peaks.add(p);
//				}
//
//				/* Save the peaks read from the file */
//				spectrumUnderConstruction.setPeaks(peaks);
//
//				spectrumUnderConstruction.setFile(inFile);
//				
//				/* clean the peaks */
//				cleanPeaks(spectrumUnderConstruction);
//
//				//Store the spectrum into output
//				out.add(spectrumUnderConstruction);

			}//if size == 2


		}//While iterator hasNext
		
		U.p("finished with mzML file: " + inFile.getName());

		return out;
	}//loadMZML Spectra

	private static ArrayList<Spectrum> loadMGFSpectra(File inFile){
		ArrayList<Spectrum> out = new ArrayList<Spectrum>();
		Spectrum temp = new Spectrum();
		
		
		//Parse the file and create a series of spectrum to add to the output
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
				
				/* begin new spectrum */
				if(line.equals("BEGIN IONS")){
					line = s.nextLine();
					
					/* header variables */
					String title = "";
					String scans = "";
					String rtInSeconds = "";
					String charge = "";
					String pepmass = "";
					
					/* peak list */
					ArrayList<Peak> peaks = new ArrayList<Peak>();
					
					while (!line.equals("END IONS")) {
						/* skip comments */
						if(line.startsWith("#")){
							line = s.nextLine();
							continue;
						}
						
						/* if the first character is a number then it is peak data */
						if (Character.isDigit(line.charAt(0))) {
							String[] chunks = line.split("\\s");
							try {
								chunks = line.split("\\s");
								line = chunks[0] + " " + chunks[1];
								Peak p = new Peak(line);
								peaks.add(p);
							} catch (Exception e) {
								//don't add it if it's bad!
							}	
						}
						
						/* else it is header data */
						else {
							String[] chunks = line.split("=");
							
							if(chunks[0].toUpperCase().startsWith("TITLE")){

								String[] brokenTitleLine = chunks[1].split("\\s");
								 title = brokenTitleLine[0].substring(brokenTitleLine[0].indexOf(':')  + 1, brokenTitleLine[0].length());
							}
							
							if(chunks[0].toUpperCase().startsWith("SCANS")){
								scans = chunks[1];
							}
							
							if(chunks[0].toUpperCase().startsWith("RTINSECONDS")){
								rtInSeconds = chunks[1];
							}
							
							if(chunks[0].toUpperCase().startsWith("CHARGE")){
								charge = chunks[1];
								if(!charge.equals("")){
									charge = charge.substring(0, charge.length() - 1);
								}
							}
							
							if(chunks[0].toUpperCase().startsWith("PEPMASS")){
								pepmass = chunks[1];
								String [] massValues = pepmass.split("\\s");
								if (massValues.length > 1) {
									pepmass = massValues[0];
								}
							}
						}
						
						line = s.nextLine();
					}
					//Loop through and load up a spectra
					
					temp.setTitle(title);
					temp.setFile(inFile);
			
					
					/* set headers */
					
					try{
						temp.setPrecursorMZ(Double.parseDouble(pepmass));
					}catch(NumberFormatException e){
						//Ignore these if they do not fit the format, and let the default value stay in the spectrum object
					}
					try{
						temp.setCharge(Integer.parseInt(charge));
					}catch(NumberFormatException e){
						//Ignore these if they do not fit the format, and let the default value stay in the spectrum object
					}
					try{
						temp.setMass(temp.getPrecursorMZ() - Definitions.HYDROGEN_MONO);
						temp.setMass(temp.getMass() * temp.getCharge());
					}catch(NumberFormatException e){
						//Ignore these if they do not fit the format, and let the default value stay in the spectrum object
					}

					try{
						temp.setScanCount(Integer.parseInt(scans));
					}catch(NumberFormatException e){
						//Ignore these if they do not fit the format, and let the default value stay in the spectrum object
					}
					try{
						temp.setRetentTime(Double.parseDouble(rtInSeconds));
					}catch(NumberFormatException e){
						//Ignore these if they do not fit the format, and let the default value stay in the spectrum object
					}
		
					
					
					//Set the peaks to this spectrum object
					temp.setPeaks(peaks);
					
					cleanPeaks(temp);
					
					temp.setTitle(inFile.getAbsolutePath().substring(inFile.getAbsolutePath().lastIndexOf('/')));
					temp.setFile(inFile);
					
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
	 * Uses jmzReader from: http://code.google.com/p/jmzreader/wiki/Welcome
	 * @param inFile
	 * @return
	 */
	private static ArrayList<Spectrum> loadMzXMLSpectra(File inFile){
		ArrayList<Spectrum> out = new ArrayList<Spectrum>();
		
		
		//Upload the MzXML fil using the jmzReader library
		try {
			JMzReader mzxmlFile = new MzXMLFile(inFile);
			
			
			Iterator<uk.ac.ebi.pride.tools.jmzreader.model.Spectrum> iter = mzxmlFile.getSpectrumIterator();
			
			
			//Iterate through the spectra information in a file and create Spectrum objects from them.
			while(iter.hasNext()){
			
				Spectrum peppySpectrum = new Spectrum();
				
				uk.ac.ebi.pride.tools.jmzreader.model.Spectrum jmzSpec = iter.next();
				
				
				//Get the precursor charge and m/z
				if(jmzSpec.getPrecursorCharge() != null){
					peppySpectrum.setCharge(jmzSpec.getPrecursorCharge());
				}
				
				if(jmzSpec.getPrecursorMZ() != null){
					peppySpectrum.setPrecursorMZ(jmzSpec.getPrecursorMZ());
				}
				
				
				

				if(jmzSpec.getPrecursorMZ() != null && jmzSpec.getPrecursorCharge() != null){
					peppySpectrum.setMass((jmzSpec.getPrecursorMZ() - Definitions.HYDROGEN_MONO) * jmzSpec.getPrecursorCharge());
				}
				

				
				ArrayList<Peak> peaks = new ArrayList<Peak>();
				Map<Double, Double> peakList = jmzSpec.getPeakList();
			
				for(Double mz: peakList.keySet()){
					Double intensity = peakList.get(mz);
					
					try {
						Peak p = new Peak( mz.doubleValue() + " " +  intensity.doubleValue());
						peaks.add(p);
					} catch (Exception e) {
						// Ignore bad peaks
					}//catch
					
				}//Iterate through the peaks and add them to the spectrum object
			
				
				//Set meta data
				peppySpectrum.setTitle(inFile.getAbsolutePath().substring(inFile.getAbsolutePath().lastIndexOf('/')));
				
				peppySpectrum.setPeaks(peaks);
				
				peppySpectrum.setFile(inFile);
				
				cleanPeaks(peppySpectrum);
				
				out.add(peppySpectrum);
				
			}//While going through peaks
			
			
		} catch (MzXMLParsingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		//Upload the needed information from the mzXML
		
		
		
		return out;
	}//loadMzXMLSpectra
	

	/**
	 * 
	 * Uses jmzreader http://code.google.com/p/jmzreader/
	 * @param inFile
	 * @return
	 */
	private static ArrayList<Spectrum> loadMS2Spectra(File inFile){
		ArrayList<Spectrum> out = new ArrayList<Spectrum>();
		
		try {
			JMzReader ms2File = new Ms2File(inFile);
			
			Iterator<uk.ac.ebi.pride.tools.jmzreader.model.Spectrum> iter = ms2File.getSpectrumIterator();
			
			
			//Iterate through the spectra information in a file and create Spectrum objects from them.
			while(iter.hasNext()){
				
				Spectrum peppySpectrum = new Spectrum();
				
				uk.ac.ebi.pride.tools.jmzreader.model.Spectrum jmzSpec = iter.next();
			
				//Get the precursor charge and m/z
				if(jmzSpec.getPrecursorCharge() != null){
					peppySpectrum.setCharge(jmzSpec.getPrecursorCharge());
				}
				
				if(jmzSpec.getPrecursorMZ() != null){
					peppySpectrum.setPrecursorMZ(jmzSpec.getPrecursorMZ());
				}
				
				
				
				//I am not sure which of these two methods are correct.  It seems that each file needs its mass calculated in a different way.
				if(jmzSpec.getPrecursorMZ() != null  && jmzSpec.getPrecursorCharge() != null ){
					peppySpectrum.setMass((jmzSpec.getPrecursorMZ()  - Definitions.HYDROGEN_MONO ) * jmzSpec.getPrecursorCharge());
				}
				

				
				ArrayList<Peak> peaks = new ArrayList<Peak>();
				Map<Double, Double> peakList = jmzSpec.getPeakList();
			
				for(Double mz: peakList.keySet()){
					Double intensity = peakList.get(mz);
					
					try {
						Peak p = new Peak( mz.doubleValue() + " " +  intensity.doubleValue());
						peaks.add(p);
					} catch (Exception e) {
						// Ignore bad peaks
					}//catch
					
				}//Iterate through the peaks and add them to the spectrum object
			
				
				//Set meta data
				peppySpectrum.setTitle(inFile.getAbsolutePath().substring(inFile.getAbsolutePath().lastIndexOf('/')));
				
				peppySpectrum.setPeaks(peaks);
				
				peppySpectrum.setFile(inFile);
				
				cleanPeaks(peppySpectrum);
				
				out.add(peppySpectrum);
			
			
			}//while
			
		} catch (JMzReaderException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return out;
		
	}//loadMS2Spectra
	
	
	
	/**
	 * loadSpectraFromFolder takes the following steps:
	 * 1) Recursively goes through folder.
	 * 2) Finds all files that end in .dta or .pkl or .txt
	 * 3) Extracts the spectra from those files
	 * 4) Adds all those spectra to one big ArrayList and returns that.
	 * @param folder  A file object that represents the folder containing spectra for this method to use
	 * @param spectra  The ArrayList to add spectra to
	 */
	private static int loadSpectraFromFolder(File folder, ArrayList<Spectrum> spectra) { 
		File [] files = folder.listFiles();
		int spectrumTotal = 0;
		
		//Loop through and grab all spectra from the folder
		for (int i = 0; i < files.length; i++) {
			if (files[i].isHidden()) continue;
			if (files[i].isDirectory()) {
				spectrumTotal += loadSpectraFromFolder(files[i], spectra);
				continue;
			}
			if (isValidFile(files[i].getName())) {
				ArrayList<Spectrum> loadedSpectra = loadSpectraFromFile(files[i]);
				if (loadedSpectra != null) {
					spectrumTotal += loadedSpectra.size();
					spectra.addAll(loadedSpectra);
				}
			}
		}

				
		return spectrumTotal;
	}//loadSpectraFromFOlder
		
	
	
	/**
	 * This is slightly more general than the method it overloads.  It creates 
	 * a spectrum ArrayList which it passes to the other loadSpectraFromFolder
	 * @param fileName
	 * @return
	 */
	private static ArrayList<Spectrum> loadSpectraFromFolder(File inFile) {
		ArrayList<Spectrum> spectra = new ArrayList<Spectrum>();
		loadSpectraFromFolder(inFile, spectra);
		return spectra;
	}//loadSpectraFromFolder
	

	/**
	 * 
	 * This method is recursively goes through folder and
	 * finds all files that end with a spectrum file extension
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
			String fileName = files[i].getName();
			if (isValidFile(fileName)) {
				spectraFiles.add(files[i]);
			}
		}
	}
	
	
	private static boolean isValidFile(String fileName){
		
		File f = new File(fileName);
		if(f.isDirectory()){
			return false;
		}
		for(int i = 0; i < fileExtensions.length; i++){
			if(fileName.toLowerCase().endsWith(fileExtensions[i])){
				
					return true;
		
			}//if
		}//for
		
		return false;
	}//isValidFile
	
	
/*End Public Static methods for loading spectra*/
	
	
	/*Methods for cleaning peaks*/
	/**
	 * cleanPeaks ensures a spectrum has a MD5, and then keeps the strongest peak in each region and sort the peaks by mass.
	 */
	private static void cleanPeaks(Spectrum spectrum) {
		
		//before we mess with the peak data, let's make sure we have the MD5
		spectrum.setMD5(spectrum.getMD5());
		
		/* also before we mess with the peak data, get the peak intensity distributions */
//		s.calculateDistributions();
	
		
		spectrum.setPeaks(getTopPeaksPerBucket(spectrum, (int) Properties.spectrumBucketSize, Properties.spectrumBucketCount));

//		s.setPeaks(keepStrongestPeakInRegions(s));
				
		sortPeaksByMass(spectrum);
		
	}//cleanPeaks
	
	

	
	/**
	 * Get top peaks per bucket attempts to get the specified number of peaks in each bucket within a spectrum.
	 * @param s Spectrum to clean
	 * @param bucketSize Size of buckets to create within the spectrum
	 * @param peaksPerBucket Number of peaks to get per bucket
	 * @return Freshly cleaned peaks from this spectrum.
	 */
	private static ArrayList<Peak> getTopPeaksPerBucket(Spectrum s, int bucketSize, int peaksPerBucket){
		ArrayList<Peak> topPeaks = new ArrayList<Peak>();
		
		
		int bucketStart = 0;
		int bucketStop = bucketStart + bucketSize;
		
		//Ensure the top scoring ION is always used for the loop stop variable
		double topScore = s.getMass();

		
		while(bucketStart < topScore){
			topPeaks.addAll(getTopPeaks(peaksPerBucket, bucketStart, bucketStop, s.getPeaks()));
			bucketStart += bucketSize;
			bucketStop += bucketSize;
		}//while
		
		
		return topPeaks;
	}//getTopPeaksPerBucket
	
	
	
	/**
	 * getTopPeaks gets the defined # of top peaks based on intensity in a specified range, and returns the masses of these top peaks.
	 * @param numTopPeaksToGet # of peaks to get in the specified range
	 * @param start inclusive lower bound of mass to search.
	 * @param stop exclusive upper bound of mass to search.
	 * @param peaks List of real peaks from a spectra to search.
	 * @return ArrayList of Double objects containing the masses of the # of top peaks specified by the input.
	 */
	private static ArrayList<Peak> getTopPeaks(int numTopPeaksToGet, int start, int stop, ArrayList<Peak> peaks){
		ArrayList<Peak> allPeaksList = new ArrayList<Peak>();
		
		for(Peak peak: peaks){
			peak.setCompareByIntensity();
			if(peak.getMass() >= start && peak.getMass() < stop){
				allPeaksList.add(peak);
			}//if
		}//for
	
		
		Collections.sort(allPeaksList);
		//Sorted by intensity 
		
		ArrayList<Peak> topPeaks = new ArrayList<Peak>();
		
		
		//Since there are not enough peaks to cull, just return all of them
		if(allPeaksList.size() < numTopPeaksToGet){
			for(Peak peak: allPeaksList){
				topPeaks.add(peak);
			}//for
			return topPeaks;
		}//if
		
		//Select just the desired number of top peaks;
		for(int i = allPeaksList.size() - numTopPeaksToGet; i < allPeaksList.size(); i++){
			topPeaks.add(allPeaksList.get(i));
		}//for
		
		return topPeaks;
	}//getTopPeaks
	

	/**
	 * keepStrongestPeakInRegions sorts the peaks by intensity, and then ensures that know two peaks are within each other by a factor of the fragmentTolerence
	 * variable from the Properties file.  
	 */
	private static ArrayList<Peak> keepStrongestPeakInRegions(Spectrum spectrum) {
		ArrayList<Peak> peaks = spectrum.getPeaks();
		sortPeaksByIntensity(spectrum);
		int start = spectrum.getPeaks().size() - 1;
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
/*End methods for cleaning peaks*/
	
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
	
}//SpectrumLoader
