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
	public static ArrayList<ProteinModification> proteinModifications = ProteinModification.getProteinModificationsFromFile(new File("resources/protein-modifications.txt"));
	
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
	
	public static double getAminoAcidWeightMono(char c) {
		if (c =='.') return 0.0;
		if (c =='A') return 71.03711;
		if (c =='C') return 103.00919;
		if (c =='D') return 115.02694;
		if (c =='E') return 129.04259;
		if (c =='F') return 147.06841;
		if (c =='G') return 57.02146;
		if (c =='H') return 137.05891;
		if (c =='I') return 113.08406;
		if (c =='K') return 128.09496;
		if (c =='L') return 113.08406;
		if (c =='M') return 131.04049;
		if (c =='N') return 114.04293;
		if (c =='P') return 97.05276;
		if (c =='Q') return 128.05858;
		if (c =='R') return 156.10111;
		if (c =='S') return 87.03203;
		if (c =='T') return 101.04768;
		if (c =='V') return 99.06841;
		if (c =='W') return 186.07931;
		if (c =='Y') return 163.06333;
		return -1;
	}
	
	public static double getAminoAcidWeightAverage(char c) {
		if (c =='.') return 0.0; 
		if (c =='A') return 71.0788; 
		if (c =='C') return 103.1448; 
		if (c =='D') return 115.0886; 
		if (c =='E') return 129.1155; 
		if (c =='F') return 147.1766; 
		if (c =='G') return 57.052; 
		if (c =='H') return 137.1412; 
		if (c =='I') return 113.1595; 
		if (c =='K') return 128.1742; 
		if (c =='L') return 113.1595; 
		if (c =='M') return 131.1986; 
		if (c =='N') return 114.1039; 
		if (c =='P') return 97.1167; 
		if (c =='Q') return 128.1308; 
		if (c =='R') return 156.1876; 
		if (c =='S') return 87.0782; 
		if (c =='T') return 101.1051; 
		if (c =='V') return 99.1326; 
		if (c =='W') return 186.2133; 
		if (c =='Y') return 163.176; 
		return -1;
	}
}
