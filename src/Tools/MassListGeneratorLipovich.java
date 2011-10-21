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
public class MassListGeneratorLipovich {
	
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
			Hashtable<Double, String> masses = new Hashtable<Double, String>();
			
			/* get the peptides for each sequence and add the calculated m/z masses to the list */
			for (Sequence sequence: sequences) {
				U.p("working on sequence " + sequence.getSequenceFile().getName());
				
				ArrayList<Peptide> peptides = sequence.extractAllPeptides(false);
				
				/* trim the peptides down to unique peptides */
				Collections.sort(peptides);
				Peptide previousPeptide = peptides.get(0);
				Peptide presentPeptide;
				for (int i = 1; i < peptides.size(); i++) {
					presentPeptide = peptides.get(i);
					if (previousPeptide.equals(presentPeptide)) {
						peptides.remove(i);
						i--;
					} else {
						previousPeptide = presentPeptide;
					}
				}
				
				String peptideString;
				double mass;
				for (Peptide peptide: peptides) {
					mass = peptide.getMass();
					
					/* ignoring peptides containing PTM-prone amino acids */
					peptideString = peptide.getAcidSequenceString();
					if (hasChar('M', peptideString) || hasChar('S', peptideString) || hasChar('T', peptideString) || hasChar('Y', peptideString)) {

					} else {
						/* adding the mass to our list */
						appendStringInHash(mass, peptide.getAcidSequenceString(), masses);
					}
				}

				peptides = null;
				System.gc();
				
				U.p(sequence.getSequenceFile().getName() + " digested.");
			}
			/* save all possible charges in one list */
			Hashtable<Double, String> allMOverZ = new Hashtable<Double, String>();
			
			/* create the files for the different charges */
			for (int charge = 1; charge <=2; charge++) {
			
				/* print each mass */
				ArrayList<Double> massList = Collections.list((masses.keys()));
				Collections.sort(massList);
				double mOverZ;
				for (Double mass: massList) {
					mOverZ = (mass / charge) + Definitions.HYDROGEN_MONO;
					String value = masses.get(mass);
					if (mOverZ < Properties.peptideMassMaximum && mOverZ > Properties.peptideMassMinimum) {

						/* add this new m/z to our big list */
						appendStringInHash(new Double(mOverZ), value + ": carge " + charge, allMOverZ);
					}
				}

			}
			
			/* save the full list of m over z */
			PrintWriter mOverZFile = new PrintWriter(new BufferedWriter(new FileWriter(new File(peptidesFolder,U.getFileNameWithoutSuffix(Properties.sequenceDirectoryOrFile) + " mOverZ.txt"))));
			PrintWriter peptideFile = new PrintWriter(new BufferedWriter(new FileWriter(new File(peptidesFolder,U.getFileNameWithoutSuffix(Properties.sequenceDirectoryOrFile) + " mOverZ and peptide.txt"))));
			
			/* print each mass */
			ArrayList<Double> mOverZList = Collections.list(allMOverZ.keys());
			Collections.sort(mOverZList);
			for (Double key: mOverZList) {
				String value = allMOverZ.get(key);
				mOverZFile.println(key);
				peptideFile.println(key + ", " + value);

			}
			mOverZFile.flush();
			mOverZFile.close();
			peptideFile.flush();
			peptideFile.close();
			
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
	
	private static void appendStringInHash(Double mass, String value, Hashtable<Double, String> masses) {
		String tag = masses.get(mass);
		if (tag != null) {
			value = tag + ", " + value;
		} 
		masses.put(mass, value);
	}

}
