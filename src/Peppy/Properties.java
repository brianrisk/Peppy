package Peppy;
import java.io.File;

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

	//properties for spectral cleaning
	public static boolean localMaximaCleaning = false;
	public static boolean highIntensityCleaning = false;
	public static int numberOfHighIntensityPeaksToRetain = 100;
	
	//when it comes to calculating theortical peptide mass, we can use mono or average
	public static boolean useMonoMass = true;
	
	//no fragments that weigh less than this will be admitted into the fragment list
	public static double peptideMassThreshold = 500.0; //daltons
	
	//number of Missed cleavages
	public static int numberOfMissedCleavages = 2;
	
	//when comparing a spectrum to a peptide, the mass may difference by as much as this amount
	public static double spectrumToPeptideMassError = 2.0; //daltons
	
	//MSMSFit
	public static double peakDifferenceThreshold = 0.5;
	public static double peakIntensityExponent = 0.33333333;
	public static double yIonDifference = 1.0;
	public static double bIonDifference = 1.0;
	
	//matches per spectrum
	public static int maximumNumberOfMatchesForASpectrum = 5;
	
	//These Files could be directories or files.
	public static File sequenceFile = new File("sequences");
	public static File spectraFile = new File("spectra");
	
	//HMM Score
	public static File HMMScoreParametersFile = new File("resources/HMMScore/ParamFiles");
	

}