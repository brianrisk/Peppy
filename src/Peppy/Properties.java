package Peppy;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Hashtable;


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
	public static int minimumNumberOfPeaksForAValidSpectrum = 10;
	
	//when it comes to calculating theoretical peptide mass, we can use mono or average
	public static boolean useMonoMass = true;
	
	//Sequence digestion
	public static int numberOfMissedCleavages = 1;
	public static double peptideMassMinimum = 500.0;
	public static double peptideMassMaximum = 10000.0;
	public static int minPeptideLength = 5;
	public static int maxPeptideLength = 80;
	public static boolean useSequenceRegion = false;
	public static int sequenceRegionStart = 0;
	public static int sequenceRegionStop = 0;
	
	//Segmenting up job for memory management
	public static int digestionWindowSize = 10000000;
	public static int desiredPeptideDatabaseSize = 10000000;
	public static int maxNumberOfProteinsToLoadAtOnce = 50000;
	
	//Splicing?
	public static boolean useSpliceVariants = false;
	
	
	//mass error tolerances in PPM
	public static double precursorTolerance = 100;
	public static double fragmentTolerance = 300;
	
	/* ion types */
//	public static double rightIonDifference = 1.0072764668; //x, y, z ion
//	public static double leftIonDifference = 1.0072764668;  //a, b, c ion
	public static double rightIonDifference = 1.0; //x, y, z ion
	public static double leftIonDifference = 1.0;  //a, b, c ion
	
	public static double minimumScore = 15;
	
	//This could be a directory or a file
	public static File sequenceDirectoryOrFile = new File("sequences");
	
	/* a list of peptide databases that will be iterated through in our search */
	public static ArrayList<File> sequenceDirectoryOrFileList;
	/* an ordered list of the database type (DNA, protein etc) for our peptide sources */
	public static ArrayList<Boolean> isSequenceFileDNAList;
	/* a list of spectra sources that will be iterated through in our search */
	public static ArrayList<File> spectraDirectoryOrFileList;
	
	//This could be a directory or a file
	public static File spectraDirectoryOrFile = new File("spectra");
	
	//FASTA files can be either DNA or amino acid sequences
	public static boolean isSequenceFileDNA = true;
	
	//Is the sequence file supposed to be read only in the first frame?
	public static boolean useOnlyForwardsFrames = false;
	
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
	

	/* for testing purposes */
	public static boolean thisIsATest = false;
	public static File testSequence;
	public static boolean testSequenceIsProtein = true;
	public static File testDirectory; 
	
	/* for custom jobs... */
	public static boolean isYale = false;
	
	/* for VCF Files*/
	public static String VCFFileString;
	
	/* FDR false discovery rate */
	public static int numberOfSpectraToUseForFDR = 10000;
	public static double maximumFDR = 0.01;
	
	/* PTMs */
	public static boolean multipass = false;
	public static int numberOfRegionsToKeep = 1000;
	public static boolean searchModifications = false;
	public static double modificationLowerBound = -0.3;
	public static double modificationUpperBound = 100;
	
	/* custom PTMs*/
	public static boolean cysteineCarbamylation = false;
	public static boolean methionineOxidation = false;
	public static boolean iodoacetamideDerivative = true;
	
	
	/* how we format our percents */
	public static NumberFormat nfPercent = NumberFormat.getPercentInstance();
	
	
	public static void loadProperties(String fileName) {
		File propertiesFile = new File(fileName);
		loadProperties(propertiesFile);
	}
	
	
	/**
	 * 
	 * @param fileName the name of our properties file
	 */
	public static void loadProperties(File propertiesFile) {
		
		/* All arrays must be cleared for multiple jobs to work!  clearing out our arrays */
		sequenceDirectoryOrFileList = new ArrayList<File>();
		isSequenceFileDNAList = new ArrayList<Boolean>();
		spectraDirectoryOrFileList = new ArrayList<File>();
		
		/* loading in the values from the properties file */
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
		
		/* for our formatting */
		nfPercent.setMaximumFractionDigits(2);
		
		
		
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
		
		/* getting the property name and the property value */
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
		if (propertyName.equals("minPeptideLength")) 
			minPeptideLength =Integer.valueOf(propertyValue);
		if (propertyName.equals("maxPeptideLength")) 
			maxPeptideLength =Integer.valueOf(propertyValue);
		
		//job parsing for memory management
		if (propertyName.equals("digestionWindowSize")) 
			digestionWindowSize =Integer.valueOf(propertyValue);
		if (propertyName.equals("desiredPeptideDatabaseSize")) 
			desiredPeptideDatabaseSize =Integer.valueOf(propertyValue);
		if (propertyName.equals("maxNumberOfProteinsToLoadAtOnce")) 
			maxNumberOfProteinsToLoadAtOnce =Integer.valueOf(propertyValue);
		
		
		if (propertyName.equals("minimumScore")) 
			minimumScore = Double.valueOf(propertyValue);
			
		
		//splicing
		if (propertyName.equals("useSpliceVariants"))
			useSpliceVariants = Boolean.valueOf(propertyValue);
		
		
		
		//spectrum cleaning		
		if (propertyName.equals("minimumNumberOfPeaksForAValidSpectrum")) 
			minimumNumberOfPeaksForAValidSpectrum =Integer.valueOf(propertyValue);
		
		
	
		if (propertyName.equals("sequenceDirectoryOrFile")) {
			sequenceDirectoryOrFile = new File(propertyValue);
			sequenceDirectoryOrFileList.add(sequenceDirectoryOrFile);
		}
		if (propertyName.equals("spectraDirectoryOrFile")) { 
			spectraDirectoryOrFile = new File(propertyValue);
			spectraDirectoryOrFileList.add(spectraDirectoryOrFile);
		}
		if (propertyName.equals("isSequenceFileDNA")) {
			isSequenceFileDNA = Boolean.valueOf(propertyValue);
			isSequenceFileDNAList.add(isSequenceFileDNA);
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
		
		/* for testing purposes */
		if (propertyName.equals("thisIsATest")) 
			thisIsATest = Boolean.valueOf(propertyValue);
		if (propertyName.equals("testSequence")) 
			testSequence = new File(propertyValue);
		if (propertyName.equals("testSequenceIsProtein")) 
			testSequenceIsProtein = Boolean.valueOf(propertyValue);
		if (propertyName.equals("testDirectory")) 
			testDirectory = new File(propertyValue);
		if (propertyName.equals("VCFFileString")) 
			VCFFileString = propertyValue;
		
		/* custom jobs */
		if (propertyName.equals("isYale")) 
			isYale = Boolean.valueOf(propertyValue);
		
		/* FDR */
		if (propertyName.equals("numberOfSpectraToUseForFDR"))
			numberOfSpectraToUseForFDR = Integer.valueOf(propertyValue);	
		if (propertyName.equals("maximumFDR"))
			maximumFDR = Double.valueOf(propertyValue);	
		
		
		
		/* multipass */
		if (propertyName.equals("multipass")) 
			multipass = Boolean.valueOf(propertyValue);
		if (propertyName.equals("numberOfRegionsToKeep"))
			numberOfRegionsToKeep = Integer.valueOf(propertyValue);	
		if (propertyName.equals("searchModifications")) 
			searchModifications = Boolean.valueOf(propertyValue);
		if (propertyName.equals("modificationLowerBound")) 
			modificationLowerBound = Double.valueOf(propertyValue);
		if (propertyName.equals("modificationUpperBound")) 
			modificationUpperBound = Double.valueOf(propertyValue);
		
		/* custom PTMs */
		if (propertyName.equals("cysteineCarbamylation")) 
			cysteineCarbamylation = Boolean.valueOf(propertyValue);
		if (propertyName.equals("methionineOxidation")) 
			methionineOxidation = Boolean.valueOf(propertyValue);
		if (propertyName.equals("iodoacetamideDerivative")) 
			iodoacetamideDerivative = Boolean.valueOf(propertyValue);
		
	}


	public static void generatePropertiesFile(File reportDir) {	
		reportDir.mkdirs();
		//set up our main index file
		File ppropertiesFile = new File(reportDir, "properties.txt");
		PrintWriter pw;
		try {
			pw = new PrintWriter(new BufferedWriter(new FileWriter(ppropertiesFile)));
			
			pw.println("##This could be a directory or a file ");
			pw.println("sequenceDirectoryOrFile " + Properties.sequenceDirectoryOrFile);
			pw.println("isSequenceFileDNA " + Properties.isSequenceFileDNA);
			pw.println("useOnlyForwardsFrames " + Properties.useOnlyForwardsFrames);
			pw.println();
			pw.println("##This could be a directory or a file ");
			pw.println("spectraDirectoryOrFile " + Properties.spectraDirectoryOrFile);
			pw.println();
			pw.println("minimumNumberOfPeaksForAValidSpectrum " + Properties.minimumNumberOfPeaksForAValidSpectrum);
			pw.println();
			pw.println("##Scoring Method ");
			pw.println("scoringMethodName " + Properties.scoringMethodName);
			pw.println();
			pw.println("##digest only part of a sequence ");
			pw.println("useSequenceRegion " + Properties.useSequenceRegion);
			pw.println("sequenceRegionStart " + Properties.sequenceRegionStart);
			pw.println("sequenceRegionStop " + Properties.sequenceRegionStop);
			pw.println();
			pw.println("##Memory contols ");
			pw.println("digestionWindowSize " + Properties.digestionWindowSize);
			pw.println("desiredPeptideDatabaseSize " + Properties.desiredPeptideDatabaseSize);
			pw.println("minimumScore " + Properties.minimumScore);
			pw.println();
			pw.println("##error thresholds in PPM");
			pw.println("precursorTolerance " + Properties.precursorTolerance);
			pw.println("fragmentTolerance " + Properties.fragmentTolerance);
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
			pw.println("##where we store our reports ");
			pw.println("reportDirectory " + Properties.reportDirectory);
			pw.println();
			pw.println("##VCF");
			pw.println("VCFFileString " + VCFFileString);
			pw.println();
			pw.println("##False Discovery Rates");
			pw.println("numberOfSpectraToUseForFDR " + numberOfSpectraToUseForFDR);
			pw.println("maximumFDR " + maximumFDR);
			pw.println();
			pw.println("##multipass");
			pw.println("multipass " + multipass);
			pw.println("searchModifications " + searchModifications);
			pw.println("modificationLowerBound " + modificationLowerBound);
			pw.println("modificationUpperBound " + modificationUpperBound);
			pw.println("numberOfRegionsToKeep " + numberOfRegionsToKeep);
			
			pw.println();
			pw.println("## static PTMs");
			pw.println("cysteineCarbamylation " + cysteineCarbamylation);
			pw.println("methionineOxidation " + methionineOxidation);
			pw.println("iodoacetamideDerivative " + iodoacetamideDerivative);
			

			
			pw.println();
			
			
	
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}	
	}
	

}