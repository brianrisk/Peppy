package Peppy;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import Math.MassError;

import uk.ac.ebi.jmzml.model.mzml.BinaryDataArray;
import uk.ac.ebi.jmzml.model.mzml.BinaryDataArrayList;
import uk.ac.ebi.jmzml.model.mzml.CVParam;
import uk.ac.ebi.jmzml.model.mzml.ParamGroup;
import uk.ac.ebi.jmzml.model.mzml.Precursor;
import uk.ac.ebi.jmzml.model.mzml.PrecursorList;
import uk.ac.ebi.jmzml.model.mzml.ScanList;
import uk.ac.ebi.jmzml.model.mzml.SelectedIonList;
import uk.ac.ebi.jmzml.xml.io.MzMLObjectIterator;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.JMzReaderException;
import uk.ac.ebi.pride.tools.ms2_parser.Ms2File;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLParsingException;


/**
 * 
 * @author Brian Risk 
 * @author David "Corvette" Thomas
 *
 */
public class SpectrumLoader {

	
	//A list of accepted file extensions peppy takes in.  When this is updated don't forget to update loadSpectra
	private static String[] fileExten = {".dta", ".pkl", ".mzml.xml", ".mzml", ".mgf", ".mz.xml", ".mzxml", ".mzxml.xml", ".ms2", ".zip"};
	
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

		ArrayList<Spectrum> spectra = new ArrayList<Spectrum>();
		
		//Load spectra based on whether it is a single file or a directory
		if (Properties.spectraDirectoryOrFile.isFile()) {
			spectra =  loadSpectra(Properties.spectraDirectoryOrFile);
		} else {
			U.p("loaded this many spectra: " + loadSpectraFromFolder(Properties.spectraDirectoryOrFile, spectra ));
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
	 * There may be more than
	 * one spectrum in a file, so this method returns an ArrayList.
	 * 
	 *This is the workhorse method for loading spectra.  In the end every spectrum is loaded using this method
	 * @param inFile
	 * @return An ArrayList of all the spectra
	 */
	public static ArrayList<Spectrum> loadSpectra(File inFile) {
		//if the spectra are in a .mzml file, then use a different method and return;
		if(inFile.getAbsolutePath().toLowerCase().endsWith(".mzml") || inFile.getAbsolutePath().toLowerCase().endsWith(".mzml.xml") ){
			return loadMZMLSpectra(inFile);
		}//if
		
		if(inFile.getAbsolutePath().toLowerCase().endsWith(".mgf")){
			return loadMGFSpectra(inFile);
		}
		
		if(inFile.getAbsolutePath().toLowerCase().endsWith(".dta")){
			return loadDTASpectra(inFile);
		}
		
		if(inFile.getAbsolutePath().toLowerCase().endsWith(".pkl")){
			return loadPKLSpectra(inFile);
		}
		if(inFile.getAbsolutePath().toLowerCase().endsWith(".mz.xml") || inFile.getAbsolutePath().toLowerCase().endsWith(".mzxml")
					|| inFile.getAbsolutePath().toLowerCase().endsWith(".mzxml.xml")){
			return loadMzXMLSpectra(inFile);
		}
		
		if(inFile.getAbsolutePath().toLowerCase().endsWith(".ms2")){
			return loadMS2Spectra(inFile);
		}
		
//		if(inFile.getAbsolutePath().toLowerCase().endsWith(".zip")){
//	
//			String tempPrefix = "tempStorage/";
//			//Name files dynamically to ensure that two zips in the same directory do not write over one another when extracted
//			String tempLoc;
//			if(inFile.getPath().contains(tempPrefix)){
//				tempLoc = inFile.getPath() +"_" + tempPrefix; //+ U.getFileNameWithoutSuffix(inFile) + "/";
//			}else{
//				tempLoc = tempPrefix + U.getFileNameWithoutSuffix(inFile) + "/";
//			}
//			
//			
//			//Uncompress and extract all files form the zip file
//			try {
//				File f = new File("ERROR");
//				File outDir = new File(tempLoc);
//				outDir.mkdirs();
//				ZipFile zipFile = new ZipFile(inFile);
//				Enumeration entries = zipFile.entries();
//				while(entries.hasMoreElements()) {
//					ZipEntry entry = (ZipEntry)entries.nextElement();
//					if(entry.isDirectory()) {
//						
//						(new File(tempLoc + entry.getName())).mkdirs();
//						continue;
//					}
//					
//					
//					String entryDirectory = "";
//					if(entry.getName().lastIndexOf('/') != -1){
//					entryDirectory = entry.getName().substring(0, entry.getName().lastIndexOf('/') + 1);
//					}
//					f = new File(tempLoc + entryDirectory);
//					f.mkdirs();
//					copyInputStream(zipFile.getInputStream(entry), new BufferedOutputStream(new FileOutputStream( tempLoc + entry.getName())));
//				}//while
//				zipFile.close();
//			
//				//Delete this zipFile if it was in the original zip file, since it will no longer need to be accessed
//				
//				//Upload all of the spectra from the newly decompressed zip file
//				ArrayList<Spectrum> out = loadSpectraFromFolder(new File(tempLoc));
//				
//
//				
//				if(inFile.getAbsolutePath().contains(tempPrefix)){
//					return out;
//				}
//
//				
//				deleteDir(new File(tempPrefix));
//				
//				//Return the results
//				return out;
//			} catch (ZipException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			} catch (IOException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//		}
		
		
		return null;
	}//loadSpectra
	
	
	/* Code to delete a directory.  Copied from: http://www.roseindia.net/tutorialhelp/comment/81393 */
	public static boolean deleteDir(File dir) {
		if (dir.isDirectory()) {
			String[] children = dir.list();
			for (int i=0; i<children.length; i++) {
				boolean success = deleteDir(new File(dir, children[i]));
				if (!success) {
					return false;
				}
			}
		}

		// The directory is now empty so delete it
		return dir.delete();
	} 
	
	public static final void copyInputStream(InputStream in, OutputStream out)
	throws IOException
	{
	byte[] buffer = new byte[1024];
	int len;
	while((len = in.read(buffer)) >= 0)
	out.write(buffer, 0, len);
	in.close();
	out.close();
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
				
				String scanStartTime = "-1";
				String scanStopTime = "-1";
				ScanList sl = spectrum.getScanList();
				if(sl != null){
					if(sl.getCvParam() != null){
						if(sl.getCvParam().iterator() != null){
							Iterator<CVParam> iter = sl.getCvParam().iterator();
							while(iter.hasNext()){
								CVParam cvp = iter.next();
								if(cvp.getUnitName() != null){
									if(cvp.getUnitName().equals("scan start time")){
										scanStartTime = cvp.getValue();
										if(cvp.getUnitName().equals("minute")){
											scanStartTime = Double.toString(Double.parseDouble(cvp.getValue()) * 60);
										}//if
									}else if(cvp.getUnitName().equals("scan stop time")){
										scanStopTime = cvp.getValue();
										if(cvp.getUnitName().equals("minute")){
											scanStopTime = Double.toString(Double.parseDouble(cvp.getValue()) * 60);
										}//if
									}//if
								}//if unName is not null
								
							}//while
						}//iterator != null
					}//cvParam != null
				}//sl != null
				
				
				temp.setScanStopTime(Double.parseDouble(scanStopTime));
				temp.setScanStartTime(Double.parseDouble(scanStartTime));
				//Get the precursor information if this spectrum has it
				//Should this be modified to require precursor information in a spectrum.  right now a spectra cannot have it and it will still be created.
				if(pcl != null){
					List<Precursor> lpc = pcl.getPrecursor();
					SelectedIonList sil = lpc.get(0).getSelectedIonList();
					List<ParamGroup> lsi = sil.getSelectedIon();
					List<CVParam> lcv = lsi.get(0).getCvParam();
					//This is the same as .PKL format files
					//0: This grabs selected ion m/z
					//1: This grabs peak intensity
					//2: This grabs charge state
					
					for(int i = 0; i < lcv.size(); i++){
						if(lcv.get(i).getName().equals("selected ion m/z")){
							temp.setPrecursorMZ(Double.parseDouble(lcv.get(i).getValue()));
						}
						if(lcv.get(i).getName().equals("charge state")){
							temp.setCharge(Integer.parseInt(lcv.get(i).getValue()));
						}
					}
					
					temp.setMass(temp.getPrecursorMZ() - Definitions.HYDROGEN_MONO);
					//Try mutliplyng by charge
					temp.setMass(temp.getMass() * temp.getCharge());
					
				}//if pcl is not null(The spectra has precursor information
			
				
				//scan and retention times
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
			  
			  temp.setFile(inFile);
			  
			  //Store the spectrum into output
			  out.add(temp);
			  
		  }//if size == 2
		  
		}//While iterator hasNext
		
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
				
				if(line.equals("BEGIN IONS")){
					line = s.nextLine();
					//Loop through and load up a spectra
					
					String[] chunks = line.split("\\s");
					
					/*Parse the header lines from this spectra*/
					String title = "";

					String scans = "";
					String RTINSECONDS = "";
					
					
					String charge = "";
					String pepmass = "";

					while(chunks.length == 1){
						//Its the title line
						chunks = line.split("=");
	
						if(chunks[0].toUpperCase().startsWith("TITLE")){

							String[] brokenTitleLine = chunks[1].split("\\s");
							 title = brokenTitleLine[0].substring(brokenTitleLine[0].indexOf(':')  + 1, brokenTitleLine[0].length());
						}
						
						if(chunks[0].toUpperCase().startsWith("SCANS")){

							scans = line.split("=")[1];
						}
						
						if(chunks[0].toUpperCase().startsWith("RTINSECONDS")){

							RTINSECONDS = line.split("=")[1];
						}
						
						if(chunks[0].toUpperCase().startsWith("CHARGE")){

							charge = line.split("=")[1];
						}
						
						if(chunks[0].toUpperCase().startsWith("PEPMASS")){
							pepmass = line.split("=")[1];
						}
						
						line = s.nextLine();
						chunks = line.split("\\s");
						
					}//while
					
					temp.setTitle(title);
					temp.setFile(inFile);
					
			
					if(!charge.equals("")){
						charge = charge.substring(0, charge.length() - 1);
					}
					//Add in precursor
					
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
						temp.setRetentTime(Double.parseDouble(RTINSECONDS));
					}catch(NumberFormatException e){
						//Ignore these if they do not fit the format, and let the default value stay in the spectrum object
					}
		
					
					
					

					
					ArrayList<Peak> peaks = new ArrayList<Peak>();
					while(s.hasNextLine()){
						
						if(line.equals("END IONS")){
							break;
						}
						
						//Only add value peaks that have a positive mass
				
						try {
							
							chunks = line.split("\\s");
							line = chunks[0] + " " + chunks[1];
							Peak p = new Peak(line);
							peaks.add(p);
						} catch (Exception e) {
							//don't add it if it's bad!
						}	
						line = s.nextLine();
						
					}//while
					
					
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
	public static int loadSpectraFromFolder(File folder, ArrayList<Spectrum> spectra) { 
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
				ArrayList<Spectrum> loadedSpectra = loadSpectra(files[i]);
				spectrumTotal += loadedSpectra.size();
				spectra.addAll(loadedSpectra);
			}
		}
		//giving all spectra an ID
		/*
		 * what the fuck is this?  why are we, inside of a recursive function, setting the ID?
		 */
		for (int i = 0; i < spectra.size(); i++) {
			spectra.get(i).setId(i);
		}
		return spectrumTotal;
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
		for(int i = 0; i < fileExten.length; i++){
			if(fileName.toLowerCase().endsWith(fileExten[i])){
				
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
	public static void cleanPeaks(Spectrum s) {
		
		//before we mess with the peak data, let's make sure we have the MD5
		s.setMD5(s.getMD5());
	
		if (Properties.maximumNumberOfPeaksforASpectrum > 0) {
			s.setPeaks(keepMostIntensePeaks(s));
		}
		
		s.setPeaks(getTopPeaksPerBucket(s, 100, 10));

		s.setPeaks(keepStrongestPeakInRegions(s));
		sortPeaksByMass(s);
		
	}//cleanPeaks
	

	
	private static ArrayList<Peak> keepMostIntensePeaks(Spectrum s) {
		ArrayList<Peak> peaks = s.getPeaks();
		if (peaks.size() <= Properties.maximumNumberOfPeaksforASpectrum) return peaks;
		sortPeaksByIntensity(s);
		while (peaks.size() > Properties.maximumNumberOfPeaksforASpectrum) {
			peaks.remove(0);
		}
		return peaks;
	}
	

	
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
		
		for(Peak p: peaks){
			p.setCompareByIntensity();
			if(p.getMass() >= start && p.getMass() < stop){
				allPeaksList.add(p);
			}//if
		}//for
	
		
		Collections.sort(allPeaksList);
		//Sorted by intensity 
		
		ArrayList<Peak> topPeaks = new ArrayList<Peak>();
		
		
		//Since there are not enough peaks to cull, just return all of them
		if(allPeaksList.size() < numTopPeaksToGet){
			for(Peak p: allPeaksList){
				topPeaks.add(p);
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
