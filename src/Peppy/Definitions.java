package Peppy;

import java.io.File;
import java.util.ArrayList;

/**
 * A class where we store all of our constants. All variables
 * should be declared "final".  All methods should be static
 * 
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class Definitions {	
	
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
	
	public final static double ITRAQ_REAGENT = 144.1021;
	
	//protein modifications
	public static ArrayList<ModificationEntry> modificationEntries = ModificationEntry.getProteinModificationsFromFile(new File("resources/protein-modifications.txt"));
	
	
}
