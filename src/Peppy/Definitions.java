package Peppy;

import java.io.File;
import java.util.ArrayList;

/**
 * A class where we store all of our constants. All variables
 * should be declared "final".  All methods should be static
 * 
 * Copyright 2012, Brian Risk
 * Released under the Netscape Public License
 * 
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
	public final static double PROTON_MONO = 1.0072764668;
	public final static double CARBON_MONO = 12.00000000;
	public final static double CARBON_AVERAGE = 12.01073590;
	public final static double WATER_MONO = OXYGEN_MONO + HYDROGEN_MONO + HYDROGEN_MONO;
	public final static double WATER_AVERAGE = OXYGEN_AVERAGE + HYDROGEN_AVERAGE + HYDROGEN_AVERAGE;
	
	//protein modifications
	public static ArrayList<ModificationEntry> modificationEntries = ModificationEntry.getProteinModificationsFromFile(new File("resources/protein-modifications.txt"));
	
	
}
