package Peppy;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import Utilities.U;

/**
 * This is were property defaults are defined.
 * Also property values are managed here.
 * If we were to import a property file, this is where the code would be.
 * @author Brian Risk
 *
 */
public class Properties {
	
	//how many processors does your computer have?  This number should be that number.
	public static int numberOfThreads = 16;
	
	//which scoring mechanism to use?
	public final static int DEFAULT_SCORE_TANDEM_FIT = 0;
	public final static int DEFAULT_SCORE_HMM = 1;
	public static int defaultScore = DEFAULT_SCORE_TANDEM_FIT;

	//properties for spectral cleaning
	public static boolean localMaximaCleaning = false;
	public static boolean highIntensityCleaning = false;
	public static int numberOfHighIntensityPeaksToRetain = 100;
	
	//when it comes to calculating theoretical peptide mass, we can use mono or average
	public static boolean useMonoMass = true;
	
	//Sequence digestion
	public static int numberOfMissedCleavages = 2;
	public static boolean onlyUsePeptidesInOpenReadingFrames = true;
	public static double peptideMassThreshold = 500.0;
	
	//when comparing a spectrum to a peptide, the mass may difference by as much as this amount
	//units are daltons.
	public static double spectrumToPeptideMassError = 2.0;
	
	//TandemFit
	public static double peakDifferenceThreshold = 0.5;
	public static double peakIntensityExponent = 0.33333333;
	public static double rightIonDifference = 1.0; //x, y, z ion
	public static double leftIonDifference = 1.0;  //a, b, c ion
	
	//matches per spectrum
	public static int maximumNumberOfMatchesForASpectrum = 5;
	
	//e value cut off
	public static double eValueCutOff = 0.01;
	
	//This could be a directory or a file
	public static File sequenceDirectoryOrFile = new File("sequences");
	
	//This could be a directory or a file
	public static File spectraDirectoryOrFile = new File("spectra");
	
	//FASTA files can be either DNA or amino acid sequences
	public static boolean isSequenceFileDNA = true;
	
	//HMM Score parameter file folder
	public static File HMMScoreParametersFile = new File("resources/HMMScore/ParamFiles");
	
	//where we store our reports
	public static File reportDirectory = new File("reports");
	
	//where we put our validation report
	public static File validationDirectory = new File("validation");
	
	public static boolean reduceDuplicateMatches = false;
	
	
	/**
	 * 
	 * @param fileName the name of our properties file
	 */
	public static void loadProperties(String fileName) {
		File file = new File(fileName);
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			String line = br.readLine();
			while (line != null) {
				setPropertyFromString(line);
				line = br.readLine();
			}
			br.close();
		} catch (FileNotFoundException e) {
			U.p("Could not find the properties file: " + fileName);
			U.p("Using default properties...");
		} catch (IOException e) {
			U.p("Could not read the properties file: " + fileName);
			e.printStackTrace();
		}
		
	}
	
	private static void setPropertyFromString(String line) {
		line = line.trim();
		if (line.equals("")) return;
		if (line.startsWith("//")) return;
		if (line.startsWith("#")) return;
		if (line.indexOf(" ") == -1) return;
		String propertyName = line.substring(0, line.indexOf(" "));
		String propertyValue = line.substring(line.indexOf(" ") + 1, line.length());
		if (propertyName.equals("numberOfThreads")) {
			numberOfThreads = Integer.valueOf(propertyValue);
		}
		
		//sequence digestion
		if (propertyName.equals("numberOfMissedCleavages")) {
			numberOfMissedCleavages =Integer.valueOf(propertyValue);
		}
		if (propertyName.equals("onlyUsePeptidesInOpenReadingFrames")) {
			onlyUsePeptidesInOpenReadingFrames = Boolean.valueOf(propertyValue);
		}
		if (propertyName.equals("peptideMassThreshold")) {
			peptideMassThreshold = Double.valueOf(propertyValue);
		}
		
		//spectrum cleaning
		if (propertyName.equals("localMaximaCleaning")) {
			localMaximaCleaning = Boolean.valueOf(propertyValue);
		}
		if (propertyName.equals("highIntensityCleaning")) {
			highIntensityCleaning = Boolean.valueOf(propertyValue);
		}
		
		if (propertyName.equals("maximumNumberOfMatchesForASpectrum")) {
			maximumNumberOfMatchesForASpectrum = Integer.valueOf(propertyValue);
		}
		if (propertyName.equals("sequenceDirectoryOrFile")) {
			sequenceDirectoryOrFile = new File(propertyValue);
		}
		if (propertyName.equals("spectraDirectoryOrFile")) {
			spectraDirectoryOrFile = new File(propertyValue);
		}
		if (propertyName.equals("isSequenceFileDNA")) {
			isSequenceFileDNA = Boolean.valueOf(propertyValue);
		}
		
		if (propertyName.equals("defaultScore")) {
			//Default to tandemFit:
			defaultScore = DEFAULT_SCORE_TANDEM_FIT;
			if (propertyValue.equals("HMM_Score")) {
				defaultScore = DEFAULT_SCORE_HMM;
			}
		}
		
		if (propertyName.equals("leftIonDifference")) {
			leftIonDifference = Double.valueOf(propertyValue);
		}
		if (propertyName.equals("rightIonDifference")) {
			rightIonDifference = Double.valueOf(propertyValue);
		}
		
		if (propertyName.equals("eValueCutOff")) {
			eValueCutOff = Double.valueOf(propertyValue);
		}
		
		
	}
	

}