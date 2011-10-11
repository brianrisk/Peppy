package Peppy;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;

import Utilities.U;

/**
 * This is were property defaults are defined.
 * Also property values are managed here.
 * If we were to import a property file, this is where the code would be.
 * @author Brian Risk
 *
 */
public class Properties {
	
	private static Hashtable<String, Property> allProperties = new Hashtable<String, Property>();
	
	//how many processors does your computer have?  This number should be that number.
	public static int numberOfThreads = Runtime.getRuntime().availableProcessors();
	
	//Define our scoring method
	public static String scoringMethodName = "Peppy.Match_IMP";
	public static MatchConstructor matchConstructor = new MatchConstructor(Properties.scoringMethodName);

	//properties for spectral cleaning
	public static boolean highIntensityCleaning = false;
	public static int numberOfHighIntensityPeaksToRetain = 100;
	public static int minimumNumberOfPeaksForAValidSpectrum = 20;
	
	//ignore spectra with large charges
	public static boolean ignoreSpectraWithChargeGreaterThanTwo = false;
	
	//when it comes to calculating theoretical peptide mass, we can use mono or average
	public static boolean useMonoMass = true;
	
	//in some experiments we will have a heavy R and K
	public static boolean useIsotopeLabeling = false;
	
	//Sequence digestion
	public static int numberOfMissedCleavages = 2;
	public static double peptideMassMinimum = 500.0;
	public static double peptideMassMaximum = 10000.0;
	public static int minPeptideLength = 5;
	public static int maxPeptideLength = 80;
	public static boolean useSequenceRegion = false;
	public static boolean useReverseDatabase = false;
	
	//Segmenting up job for memory management
	public static int numberOfSpectraPerSegment = 60000;
	public static int digestionWindowSize = 10000000;
	public static int desiredPeptideDatabaseSize = 10000000;
	public static int maxNumberOfProteinsToLoadAtOnce = 50000;
	
	//Splicing?
	public static boolean useSpliceVariants = false;
	public static int sequenceRegionStart = 0;
	public static int sequenceRegionStop = 0;
	
	//mass error tolerances
	public static double precursorTolerance = 2.0;
	public static double fragmentTolerance = 0.3;
	
	//TandemFit
	public static double peakIntensityExponent = 0.33333333;
	public static double rightIonDifference = 1.0; //x, y, z ion
	public static double leftIonDifference = 1.0;  //a, b, c ion
	public static double YBtrue = 1.1;
	public static double YBfalse = 0.9;
	public static double BYtrue = 1.1;
	public static double BYfalse = 0.9;
	
	
	//matches per spectrum
	public static int maximumNumberOfMatchesForASpectrum = 5;
	
	//cut offs.  sexy, sexy cutoffs.
	public static double maxEValue = 0.1;
	public static double maxIMP = 0.000001;
	
	//This could be a directory or a file
	public static File sequenceDirectoryOrFile = new File("sequences");
	
	//This could be a directory or a file
	public static File spectraDirectoryOrFile = new File("spectra");
	
	//FASTA files can be either DNA or amino acid sequences
	public static boolean isSequenceFileDNA = true;
	
	//do the sequence files contain many sequences?
	
	//Is the sequence file supposed to be read only in the first frame?
	public static boolean useOnlyForwardsFrames = false;
	
	//HMM Score parameter file folder
	public static File HMMScoreParametersFile = new File("resources/HMMScore/ParamFiles");
	
	//where we store our reports
	public static File reportDirectory = new File("reports");
	
	//where we put our validation report
	public static File validationDirectory = new File("validation");
	
	//Report related
	public static boolean createHTMLReport = true;
	
	public static String reportWebSuffix = ".html";
	public static File reportWebHeaderFile = new File("resources/reports/header.txt");
	public static File reportWebHeaderSubFile = new File("resources/reports/header-sub.txt");
	public static File reportWebFooterFile = new File("resources/reports/footer.txt");
	public static File reportWebTableHeader = new File("resources/reports/index-table-header.txt");
	
	//the number of nucleotides away from a specific location on a chromosome for it to be
	//considered part of the "neighborhood"
	public static int locusNeighborhood = 3000;
	
	//Multiple Modifications
	public static double multiModPrecursorMargin = 0.1;
	
	/* for testing purposes */
	public static File testSequence;
	public static boolean testSequenceIsProtein = true;
	public static File testDirectory; 
	
	/* for custom jobs... */
	public static boolean isYale = false;
	
	/*for VCF Files*/
	public static String VCFFileString;
	
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
	
	public static void setAllProperties() {
		numberOfThreads = allProperties.get("numberOfThreads").getInt();
	}
	
	private static void setPropertyFromString(String line) {
		line = line.trim();
		
		/* ignore blank lines */
		if (line.equals("")) return;
		
		/* ignore comments */
		if (line.startsWith("//")) return;
		if (line.startsWith("#")) return;
		
		/* ignore lines that do not have a space in them */
		if (line.indexOf(" ") == -1) return;
		
		/* getting the property name and the propert value */
		String propertyName = line.substring(0, line.indexOf(" "));
		String propertyValue = line.substring(line.indexOf(" ") + 1, line.length());
		
		//sequence digestion
		if (propertyName.equals("numberOfMissedCleavages")) 
			numberOfMissedCleavages =Integer.valueOf(propertyValue);
		if (propertyName.equals("peptideMassThreshold")) 
			peptideMassMinimum = Double.valueOf(propertyValue);
		if (propertyName.equals("peptideMassMinimum")) 
			peptideMassMinimum = Double.valueOf(propertyValue);
		if (propertyName.equals("peptideMassMaximum")) 
			peptideMassMaximum = Double.valueOf(propertyValue);
		if (propertyName.equals("useSequenceRegion"))
			useSequenceRegion = Boolean.valueOf(propertyValue);
		if (propertyName.equals("sequenceRegionStart")) 
			sequenceRegionStart =Integer.valueOf(propertyValue);
		if (propertyName.equals("sequenceRegionStop")) 
			sequenceRegionStop =Integer.valueOf(propertyValue);
		if (propertyName.equals("useIsotopeLabeling"))
			useIsotopeLabeling = Boolean.valueOf(propertyValue);
		if (propertyName.equals("minPeptideLength")) 
			minPeptideLength =Integer.valueOf(propertyValue);
		if (propertyName.equals("maxPeptideLength")) 
			maxPeptideLength =Integer.valueOf(propertyValue);
		if (propertyName.equals("useReverseDatabase"))
			useReverseDatabase = Boolean.valueOf(propertyValue);
		
		//job parsing for memory management
		if (propertyName.equals("numberOfSpectraPerSegment")) 
			numberOfSpectraPerSegment =Integer.valueOf(propertyValue);
		if (propertyName.equals("digestionWindowSize")) 
			digestionWindowSize =Integer.valueOf(propertyValue);
		if (propertyName.equals("desiredPeptideDatabaseSize")) 
			desiredPeptideDatabaseSize =Integer.valueOf(propertyValue);
		if (propertyName.equals("maxNumberOfProteinsToLoadAtOnce")) 
			maxNumberOfProteinsToLoadAtOnce =Integer.valueOf(propertyValue);
		
		
		//splicing
		if (propertyName.equals("useSpliceVariants"))
			useSpliceVariants = Boolean.valueOf(propertyValue);
		
		
		
		//spectrum cleaning
		if (propertyName.equals("highIntensityCleaning"))
			highIntensityCleaning = Boolean.valueOf(propertyValue);
		
		if (propertyName.equals("ignoreSpectraWithChargeGreaterThanTwo"))
			ignoreSpectraWithChargeGreaterThanTwo = Boolean.valueOf(propertyValue);
		
		if (propertyName.equals("minimumNumberOfPeaksForAValidSpectrum")) 
			minimumNumberOfPeaksForAValidSpectrum =Integer.valueOf(propertyValue);
		
		
		
		if (propertyName.equals("maximumNumberOfMatchesForASpectrum"))
			maximumNumberOfMatchesForASpectrum = Integer.valueOf(propertyValue);
		if (propertyName.equals("sequenceDirectoryOrFile"))
			sequenceDirectoryOrFile = new File(propertyValue);
		if (propertyName.equals("spectraDirectoryOrFile")) 
			spectraDirectoryOrFile = new File(propertyValue);
		if (propertyName.equals("isSequenceFileDNA")) {
			isSequenceFileDNA = Boolean.valueOf(propertyValue);
		}
		if (propertyName.equals("useOnlyForwardsFrames")) {
			useOnlyForwardsFrames = Boolean.valueOf(propertyValue);
		}
		
		
		
		
		//Scoring method
		if (propertyName.equals("scoringMethodName")) {
			scoringMethodName = propertyValue;
			matchConstructor = new MatchConstructor(scoringMethodName);
		}
		
		if (propertyName.equals("leftIonDifference"))
			leftIonDifference = Double.valueOf(propertyValue);
		if (propertyName.equals("rightIonDifference")) 
			rightIonDifference = Double.valueOf(propertyValue);
		
		//e value
		if (propertyName.equals("maxEValue")) 
			maxEValue = Double.valueOf(propertyValue);
		
		//imp
		if (propertyName.equals("maxIMP")) 
			maxIMP = Double.valueOf(propertyValue);

		//Ion and precursor threshold values
		if (propertyName.equals("precursorTolerance")) 
			precursorTolerance = Double.valueOf(propertyValue);
		if (propertyName.equals("fragmentTolerance")) 
			fragmentTolerance = Double.valueOf(propertyValue);
		
		//Old ion/precursor value names
		if (propertyName.equals("spectrumToPeptideMassError")) 
			precursorTolerance = Double.valueOf(propertyValue);
		if (propertyName.equals("peakDifferenceThreshold")) 
			fragmentTolerance = Double.valueOf(propertyValue);
		
		
		//reports
		if (propertyName.equals("reportDirectory")) 
			reportDirectory = new File(propertyValue);
		if (propertyName.equals("createHTMLReport")) 
			createHTMLReport = Boolean.valueOf(propertyValue);
		
		//Multi PTMS
		if (propertyName.equals("multiModPrecursorMargin")) 
			multiModPrecursorMargin = Double.valueOf(propertyValue);
		
		/* for testing purposes */
		if (propertyName.equals("testSequence")) 
			testSequence = new File(propertyValue);
		if (propertyName.equals("testSequenceIsProtein")) 
			testSequenceIsProtein = Boolean.valueOf(propertyValue);
		if (propertyName.equals("testDirectory")) 
			testDirectory = new File(propertyValue);
		if (propertyName.equals("VCFFileString")) 
			VCFFileString = propertyValue;
		
	}


	public static void generatePropertiesFile(File reportDir) {	
		reportDir.mkdirs();
		//set up our main index file
		File ppropertiesFile = new File(reportDir, reportDir.getName() + "_properties.txt");
		PrintWriter pw;
		try {
			pw = new PrintWriter(new BufferedWriter(new FileWriter(ppropertiesFile)));
			
			pw.println("##Properties which define how a sequence is digested");
			pw.println("isSequenceFileDNA " + Properties.isSequenceFileDNA);
			pw.println("useOnlyForwardsFrames " + Properties.useOnlyForwardsFrames);
			pw.println("useIsotopeLabeling " + Properties.useIsotopeLabeling);
			pw.println("useReverseDatabase " + Properties.useReverseDatabase);
			pw.println();
			pw.println("##This could be a directory or a file ");
			pw.println("sequenceDirectoryOrFile " + Properties.sequenceDirectoryOrFile);
			pw.println();
			pw.println("##This could be a directory or a file ");
			pw.println("spectraDirectoryOrFile " + Properties.spectraDirectoryOrFile);
			pw.println();
			pw.println("##Scoring Method ");
			pw.println("scoringMethodName " + Properties.scoringMethodName);
			pw.println();
			pw.println("##retain 100 most intense peaks");
			pw.println("highIntensityCleaning " + Properties.highIntensityCleaning);
			pw.println();
			pw.println("##digest only part of a sequence ");
			pw.println("useSequenceRegion " + Properties.useSequenceRegion);
			pw.println("sequenceRegionStart " + Properties.sequenceRegionStart);
			pw.println("sequenceRegionStop " + Properties.sequenceRegionStop);
			pw.println();
			pw.println("##limit returned matches by confidence ");
			pw.println("maxEValue " + Properties.maxEValue);
			pw.println("maxIMP " + Properties.maxIMP);
			pw.println();
			pw.println("##a preference for digestion of large DNA windows ");
			pw.println("digestionWindowSize " + Properties.digestionWindowSize);
			pw.println();
			pw.println("##number of spectra to process at once ");
			pw.println("numberOfSpectraPerSegment " + Properties.numberOfSpectraPerSegment);
			pw.println();
			pw.println("##how much precursor mass / theoretical mass difference should we tolerate? ");
			pw.println("spectrumToPeptideMassError " + Properties.precursorTolerance);
			pw.println();
			pw.println("##TandemFit property ");
			pw.println("peakDifferenceThreshold " + Properties.fragmentTolerance);
			pw.println();
			pw.println("##Report variables ");
			pw.println("createHTMLReport " + Properties.createHTMLReport);
			pw.println();
			pw.println("##no fragments that weigh less than this will be admitted into the fragment list ");
			pw.println("##units are daltons. ");
			pw.println("peptideMassMinimum " + Properties.peptideMassMinimum);
			pw.println("peptideMassMaximum " + Properties.peptideMassMaximum);
			pw.println();
			pw.println("numberOfMissedCleavages " + Properties.numberOfMissedCleavages);
			pw.println();
			pw.println("##splicing ");
			pw.println("useSpliceVariants " + Properties.useSpliceVariants);
			pw.println();
			pw.println("##This is per sequence file, so if this value is 5 and you use 7 FASTA files ");
			pw.println("##it will produce (at least) 35 matches per spectrum ");
			pw.println("##the final number of results also varies depending on the digestionWindowSize ");
			pw.println("maximumNumberOfMatchesForASpectrum " + Properties.maximumNumberOfMatchesForASpectrum);
			pw.println();
			pw.println("##where we store our reports ");
			pw.println("reportDirectory " + Properties.reportDirectory);
			pw.println();
			pw.println("##VCF");
			pw.println("VCFFileString " + VCFFileString);
			pw.println();
			pw.println();
	
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}	
	}
	

}