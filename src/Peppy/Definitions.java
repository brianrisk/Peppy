package Peppy;

import java.io.File;
import java.util.ArrayList;

/**
 * A class where we store all of our constants. All variables
 * should be declared "final".  All methods should be static
 * @author Brian Risk
 *
 */
public class Definitions {	
	public final static char[] aminoAcidList = {'K','N','K','N','T','T','T','T','R','S','R','S','I','I','M','I','Q','H','Q','H','P','P','P','P','R','R','R','R','L','L','L','L','E','D','E','D','A','A','A','A','G','G','G','G','V','V','V','V','.','Y','.','Y','S','S','S','S','.','C','W','C','L','F','L','F'};
	public final static double OXYGEN_MONO = 15.99491463;
	public final static double OXYGEN_AVERAGE = 15.99940494;
	public final static double NITROGEN_MONO = 14.00307400;
	public final static double NITROGEN_AVERAGE = 14.00674309;
	public final static double HYDROGEN_MONO = 1.00782504;
	public final static double HYDROGEN_AVERAGE = 1.00794076;
	public final static double CARBON_MONO = 12.00000000;
	public final static double CARBON_AVERAGE = 12.01073590;
	//TODO this .98 is from my observations of theoretical and measured peptide mass differences.  check validity.
	public final static double WATER_MONO = OXYGEN_MONO + HYDROGEN_MONO + HYDROGEN_MONO;
	public final static double WATER_AVERAGE = OXYGEN_AVERAGE + HYDROGEN_AVERAGE + HYDROGEN_AVERAGE;
	
	//indicies in our amino acid list which define trypsin cleavages
	public final static int [] NoCleavageBefore = {20, 23, 22, 21}; 
	public final static int [] Cleavages = {8, 10, 24, 25, 26, 0, 2}; 
	
	//protein modifications
	public static ArrayList<Modification> modifications = Modification.getProteinModificationsFromFile(new File("resources/protein-modifications.txt"));
	
	//HMM state constant
	public final static int NUMBER_OF_IONS =12;
	public final static int B_ION =0;
	public final static int Y_ION= 1;
	public final static int A_ION =2;
	public final static int J_ION = 3;
	public final static int B17_ION =4;
	public final static int Y17_ION =5;
	public final static int B18_ION =6;
	public final static int Y18_ION =7;
	public final static int IMM_ION =8;
	public final static int INTERNAL_ION = 9;
	public final static int INTERNAL17 = 10;
	public final static int INTERNAL18 =11;
	
	
}
