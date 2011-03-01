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
	public final static int DEFAULT_SCORE_IMP = 2;
	public static int defaultScore = DEFAULT_SCORE_IMP;

	//properties for spectral cleaning
	public static boolean highIntensityCleaning = false;
	public static int numberOfHighIntensityPeaksToRetain = 100;
	
	//when it comes to calculating theoretical peptide mass, we can use mono or average
	public static boolean useMonoMass = true;
	
	//Sequence digestion
	public static int numberOfMissedCleavages = 2;
	public static boolean onlyUsePeptidesInOpenReadingFrames = true;
	public static double peptideMassThreshold = 500.0;
	public static int peptideMaximumLength = 80;
	public static int digestionWindowSize = 25000000;
	
	//Splicing?
	public static boolean useSpliceVariants = false;
	public static boolean useSequenceRegion = false;
	public static int sequenceRegionStart = 0;
	public static int sequenceRegionStop = 0;
	
	//when comparing a spectrum to a peptide, the mass may difference by as much as this amount
	//units are daltons.
	public static double spectrumToPeptideMassError = 2.0;
	
	//TandemFit
	public static double peakDifferenceThreshold = 0.5;
//	public static double peakIntensityExponent = 0.33333333;
	public static double peakIntensityExponent = 0.30;
	public static double rightIonDifference = 1.0; //x, y, z ion
	public static double leftIonDifference = 1.0;  //a, b, c ion
//	public static double YBtrue = 1.1;
//	public static double YBfalse = 1.2;
//	public static double BYtrue = 1.3;
//	public static double BYfalse = 0.9;
	public static double YBtrue = 1.17;
	public static double YBfalse = 1.43;
	public static double BYtrue = 1.1;
	public static double BYfalse = 0.9;
	
	
	//matches per spectrum
	public static int maximumNumberOfMatchesForASpectrum = 5;
	
	//e value cut off
	public static double eValueCutOff = 0.03359587957603186;
	public static boolean useEValueCutOff = true;
	
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
	public static String reportDirectoryTitle = "report";
	
	//where we put our validation report
	public static File validationDirectory = new File("validation");
	
	public static boolean reduceDuplicateMatches = false;
	
	//Report related
	public static boolean createHTMLReport = true;
	public static boolean generateNeighborhoodReport = false;
	public static boolean generateSequenceReport = false;
	public static boolean generateSpectrumReport = true;
	
	public static String reportWebSuffix = ".html";
	public static File reportWebHeaderFile = new File("resources/reports/header.txt");
	public static File reportWebHeaderSubFile = new File("resources/reports/header-sub.txt");
	public static File reportWebFooterFile = new File("resources/reports/footer.txt");
	public static File reportWebTableHeader = new File("resources/reports/index-table-header.txt");
	
	//the number of nucleotides away from a specific location on a chromosome for it to be
	//considered part of the "neighborhood"
	public static int locusNeighborhood = 3000;
	
	//Multiple Modifications
	public static double multiModPrecursorMargin = 0.006;
	
	public static void loadProperties(String fileName) {
		File propertiesFile = new File(fileName);
		loadProperties(propertiesFile);
	}
	
	
	/**
	 * 
	 * @param fileName the name of our properties file
	 */
	public static void loadProperties(File propertiesFile) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(propertiesFile));
			String line = br.readLine();
			while (line != null) {
				setPropertyFromString(line);
				line = br.readLine();
			}
			br.close();
		} catch (FileNotFoundException e) {
			U.p("Could not find the properties file: " + propertiesFile.getName());
			U.p("Using default properties...");
		} catch (IOException e) {
			U.p("Could not read the properties file: " + propertiesFile.getName());
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
		if (propertyName.equals("numberOfMissedCleavages")) 
			numberOfMissedCleavages =Integer.valueOf(propertyValue);
		if (propertyName.equals("onlyUsePeptidesInOpenReadingFrames")) 
			onlyUsePeptidesInOpenReadingFrames = Boolean.valueOf(propertyValue);
		if (propertyName.equals("peptideMassThreshold")) 
			peptideMassThreshold = Double.valueOf(propertyValue);
		if (propertyName.equals("digestionWindowSize")) 
			digestionWindowSize =Integer.valueOf(propertyValue);
		
		//splicing
		if (propertyName.equals("useSpliceVariants"))
			useSpliceVariants = Boolean.valueOf(propertyValue);
		if (propertyName.equals("useSequenceRegion"))
			useSequenceRegion = Boolean.valueOf(propertyValue);
		if (propertyName.equals("sequenceRegionStart")) 
			sequenceRegionStart =Integer.valueOf(propertyValue);
		if (propertyName.equals("sequenceRegionStop")) 
			sequenceRegionStop =Integer.valueOf(propertyValue);
		
		
		//spectrum cleaning
		if (propertyName.equals("highIntensityCleaning"))
			highIntensityCleaning = Boolean.valueOf(propertyValue);
		
		
		if (propertyName.equals("maximumNumberOfMatchesForASpectrum"))
			maximumNumberOfMatchesForASpectrum = Integer.valueOf(propertyValue);
		if (propertyName.equals("sequenceDirectoryOrFile"))
			sequenceDirectoryOrFile = new File(propertyValue);
		if (propertyName.equals("spectraDirectoryOrFile")) 
			spectraDirectoryOrFile = new File(propertyValue);
		if (propertyName.equals("isSequenceFileDNA")) {
			isSequenceFileDNA = Boolean.valueOf(propertyValue);
		}
		
		if (propertyName.equals("defaultScore")) {
			if (propertyValue.equals("HMM_Score")) {
				defaultScore = DEFAULT_SCORE_HMM;
			}
			if (propertyValue.equals("IMP")) {
				defaultScore = DEFAULT_SCORE_IMP;
			}
			if (propertyValue.equals("TandemFit")) {
				defaultScore = DEFAULT_SCORE_TANDEM_FIT;
			}
		}
		
		if (propertyName.equals("leftIonDifference"))
			leftIonDifference = Double.valueOf(propertyValue);
		if (propertyName.equals("rightIonDifference")) 
			rightIonDifference = Double.valueOf(propertyValue);
		
		//e value
		if (propertyName.equals("eValueCutOff")) 
			eValueCutOff = Double.valueOf(propertyValue);
		if (propertyName.equals("useEValueCutOff")) 
			useEValueCutOff = Boolean.valueOf(propertyValue);
		
		//matches
		if (propertyName.equals("spectrumToPeptideMassError")) 
			spectrumToPeptideMassError = Double.valueOf(propertyValue);
		if (propertyName.equals("peakDifferenceThreshold")) 
			peakDifferenceThreshold = Double.valueOf(propertyValue);
		
		//reports
		if (propertyName.equals("reportDirectoryTitle")) 
			reportDirectoryTitle =propertyValue;
		if (propertyName.equals("reportDirectory")) 
			reportDirectory = new File(propertyValue);
		if (propertyName.equals("createHTMLReport")) 
			createHTMLReport = Boolean.valueOf(propertyValue);
		if (propertyName.equals("generateNeighborhoodReport")) 
			generateNeighborhoodReport = Boolean.valueOf(propertyValue);
		if (propertyName.equals("generateSequenceReport")) 
			generateSequenceReport = Boolean.valueOf(propertyValue);
		if (propertyName.equals("generateSpectrumReport")) 
			generateSpectrumReport = Boolean.valueOf(propertyValue);
		
		//Multi PTMS
		if (propertyName.equals("multiModPrecursorMargin")) 
			multiModPrecursorMargin = Double.valueOf(propertyValue);
		
		
	}
	

}