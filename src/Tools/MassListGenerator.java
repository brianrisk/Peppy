package Tools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

import Peppy.AminoAcids;
import Peppy.Definitions;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.Sequence;
import Utilities.U;

/**
 * Takes the protein file specified in your properties.txt file.
 * Digests it into peptides.
 * Makes a list of masses of those peptides.
 * Goes through the peptide list and adds to the mass list if the
 * peptide could have an oxidation or phosphorylation, where it errs
 * on the side of oxidation if both are possible.
 * Creates two more lists for charge 2 and charge 3.
 * @author Brian Risk
 *
 */
public class MassListGenerator {
	
	private static final double oxidation = 15.994915;
	private static final double phosphorylation = 79.966331;
	
	public static void main(String args[]) {
		init("properties.txt");
		
		U.p("Generating mass list");
		/* Get references to our sequence files */
		ArrayList<Sequence> sequences = Sequence.loadSequenceFiles(Properties.sequenceDirectoryOrFile);
		
		try {
			/* set up the output folder */
			File peptidesFolder = new File ("exported mass lists");
			peptidesFolder.mkdir();
			
			/* set up our hashtable for our masses */
			Hashtable<Double, Double> masses = new Hashtable<Double, Double>();
			
			/* get the peptides for each sequence and add the calculated m/z masses to the list */
			for (Sequence sequence: sequences) {
				U.p("working on sequence " + sequence.getSequenceFile().getName());
				
				ArrayList<Peptide> peptides = sequence.extractAllPeptides(false);
				String peptideString;
				double mass;
				for (Peptide peptide: peptides) {
					mass = peptide.getMass();
					
					/* adding the mass to our list */
					masses.put(mass, mass);
					
					/* calculating mass of potential PTM */
					/* if oxidation is present, we ignore the possibility of phosphorylation */
					peptideString = peptide.getAcidSequenceString();
					if (hasChar('M', peptideString)) {
						mass += oxidation;
						masses.put(mass, mass);
					} else {
						if (hasChar('S', peptideString) || hasChar('T', peptideString) || hasChar('Y', peptideString)) {
							mass += phosphorylation;
							masses.put(mass, mass);
						}
					}
				}

				peptides = null;
				System.gc();
				
				U.p(sequence.getSequenceFile().getName() + " digested.");
			}
			/* save all possible charges in one list */
			Hashtable<Double, Double> allMOverZ = new Hashtable<Double, Double>();
			
			/* create the files for the different charges */
			for (int charge = 1; charge <=4; charge++) {
				/* set up our mass list file */
				PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(peptidesFolder,U.getFileNameWithoutSuffix(Properties.sequenceDirectoryOrFile) + " charge=" + charge + ".txt"))));
				
				/* print each mass */
				ArrayList<Double> massList = new ArrayList<Double>(masses.values());
				Collections.sort(massList);
				double mOverZ;
				for (Double d: massList) {
					mOverZ = (d / charge) + Definitions.HYDROGEN_MONO;
					if (mOverZ < 2000 && mOverZ > 400) {
						pw.println(mOverZ);
						
						/* add this new m/z to our big list */
						allMOverZ.put(mOverZ, mOverZ);
					}
				}
				
				pw.flush();
				pw.close();
			}
			
			/* save the full list of m over z */
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(peptidesFolder,U.getFileNameWithoutSuffix(Properties.sequenceDirectoryOrFile) + " all charges.txt"))));
			
			/* print each mass */
			ArrayList<Double> mOverZList = new ArrayList<Double>(allMOverZ.values());
			Collections.sort(mOverZList);
			for (Double d: mOverZList) {
				pw.println(d);

			}
			pw.flush();
			pw.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		U.p("done");
	}
	
	private static boolean hasChar(char needle, String haystack) {
		return (haystack.indexOf(needle) != -1);
	}
	
	
	public static void init(String propertiesFile) {
		System.setProperty("java.awt.headless", "true"); 
		Properties.loadProperties(propertiesFile);
		AminoAcids.init();
	}

}
