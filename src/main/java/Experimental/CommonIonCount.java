package Experimental;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;

import Peppy.AminoAcids;
import Peppy.Peppy;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.Sequence;
import Peppy.U;

public class CommonIonCount {
	
	public static int rows;
	public static int cols;
	public static int[][] results;
	
	
	public static void main(String args[]) {
		U.p("beginning: CommonIonCount");
		U.startStopwatch();
		
		/* set up initial state */
		Peppy.init(args);
		
		rows = (int) (Properties.peptideMassMaximum - Properties.peptideMassMinimum + 1);
		cols = (int) ((Properties.peptideMassMaximum - Properties.peptideMassMinimum) / AminoAcids.getWeightMono(AminoAcids.G) );
		results = new int[rows][cols];
		
		/* load protein/genome sequence files */
		ArrayList<Sequence> sequences = Sequence.loadSequenceFiles(Properties.sequenceDirectoryOrFile);	
		
		collectResults(sequences);
		
		writeResults();
		
		U.stopStopwatch();
		U.p("done");
	}
	
	public static void writeResults() {
		/* write the results */
		try {
			PrintWriter pw = new PrintWriter(new FileWriter("CommonIonCount.txt"));
			for (int rowIndex = 0; rowIndex < rows; rowIndex++) {
				for (int colIndex = 0; colIndex < cols; colIndex++) {
					pw.print(results[rowIndex][colIndex]);
					if (colIndex != cols) {
						pw.print("\t");
					}
				}
				pw.println();
			}
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	private static void collectResults(ArrayList<Sequence> sequences) {
			
		/* track which sequence we are examining */
		int sequenceIndex = 0;
			
		/* loops until we have gone through all of our sequences */
		while (true) {
			
			/* Extract a decent size of peptides.  Sequences may be short, so this
			 * goes through each sequence and extracts peptides until a desired
			 * threshold has been been reached.
			 */
			
			/*  Initialize our list of peptides */
			ArrayList<Peptide> peptides = new ArrayList<Peptide>(Properties.desiredPeptideDatabaseSize);
			peptides.ensureCapacity(Properties.desiredPeptideDatabaseSize);
			
			/* This is where we get a chunk of peptides */
			ArrayList<Peptide> peptideSegment = new ArrayList<Peptide>();
			
			/* if we use a region, extract that region, else go through all sequences getting a chunk at a time. */
			if (Properties.useSequenceRegion) {
				
				/* collect peptides */
				peptideSegment = sequences.get(sequenceIndex).extractMorePeptides(false);
				peptides.addAll(peptideSegment);
				
				/* reset the sequence so that the next batch of spectra can scan it */
				sequences.get(sequenceIndex).reset();
				
			} else {
				while (peptides.size() < Properties.desiredPeptideDatabaseSize) {
					
					/* clear previous chunk of peptides and reclaim memory */
					if (peptideSegment != null) {	
						peptideSegment.clear();
						System.gc();
					}	
					
					/* collect peptides */
					peptideSegment = sequences.get(sequenceIndex).extractMorePeptides(false);
						
					/* advance to the next sequence if we don't have any more peptides in this sequence */
					if (peptideSegment == null) {
						sequences.get(sequenceIndex).reset();
						sequenceIndex++;
					
					/* add peptides to the main list if we have some to add */
					} else {
						peptides.addAll(peptideSegment);
					}
					
					/* if we have advanced past the last sequence, then exit this loop */
					if (sequenceIndex == sequences.size()) {
						break;
					}
					
				}
			}
			
			
			/* checks to see if we have any peptides in this chunk */
			if (peptides.size() > 0) {
				
				/* remove redundant peptides */
				Hashtable<String, Peptide> reducedPeptides = new Hashtable<String, Peptide>();
				for (Peptide peptide: peptides) reducedPeptides.put(peptide.getAcidSequenceString(), peptide);
				peptides = new ArrayList<Peptide>(reducedPeptides.values());
				
				/*report */
				U.p("working on peptide group totaling: " + peptides.size());
				
				comparePeptidesInArray(peptides);
				
				/* incrementally save results */
				writeResults();
				
				/* free up memory */
				peptides.clear();
				System.gc();
			
			}
				
			/* break if we have covered our last sequence or if we are only using a region */
			if (sequenceIndex == sequences.size() || Properties.useSequenceRegion) {
				break;
			}

		}
	}
	
	public static void comparePeptidesInArray(ArrayList<Peptide> peptides) {
		/* create hashtable of arrays grouped by rounded mass */
		Hashtable<Integer, ArrayList<Peptide>> peptideArrays = new Hashtable<Integer, ArrayList<Peptide>>();
		for (Peptide peptide: peptides) {
			int roundedMass = (int) Math.round(peptide.getMass());
			ArrayList<Peptide> onterPeptidesOfEqualMass = peptideArrays.get(roundedMass);
			if (onterPeptidesOfEqualMass == null) onterPeptidesOfEqualMass = new ArrayList<Peptide>();
			onterPeptidesOfEqualMass.add(peptide);
			peptideArrays.put(roundedMass, onterPeptidesOfEqualMass);
		}
		
		/* go through each of the arrays and compare all peptides within each array */
		ArrayList<ArrayList<Peptide>> collectedArrays = new ArrayList<ArrayList<Peptide>>(peptideArrays.values());
		
		IonComparisonServer ionComparisonServer = new IonComparisonServer(collectedArrays);
		ionComparisonServer.findMatches();
	}
	
	public static int compare(int peptideMassIndex, Hashtable<Integer, Integer> ionHashA, Hashtable<Integer, Integer> ionHashB) {
		
		int numberOfPeaksInCommon = 0;
		
		/* find which B and Y ions for peptideB match */
		Integer peak;
		for (Enumeration<Integer> e = ionHashA.elements() ; e.hasMoreElements() ;) {
	         peak = e.nextElement();
	         if (ionHashB.get(peak) != null) numberOfPeaksInCommon++;
	     }
		
		results[peptideMassIndex][numberOfPeaksInCommon]++;
		
		return numberOfPeaksInCommon;

	}
	
	public static Hashtable<Integer, Integer> getIons(Peptide peptide) {
		/* create a hashtable of the ions for peptideA */
//		ArrayList<Integer> peakAYIons = getYIons(peptide);
		ArrayList<Integer> peakABIons = getBIons(peptide);
		Hashtable<Integer, Integer> ionHash = new Hashtable<Integer, Integer>();
//		for (Integer peak: peakAYIons) ionHash.put(peak, peak);
		for (Integer peak: peakABIons) ionHash.put(peak, peak);
		return ionHash;
	}
	
	public static ArrayList<Integer> getYIons(Peptide peptide) {
		ArrayList<Integer> out = new ArrayList<Integer>();
		double theoreticalPeakMass = peptide.getMass() + Properties.rightIonDifference;
		byte [] acidSequence = peptide.getAcidSequence();
		for (int i = 0; i < peptide.getLengthMinusOne(); i++) {
			theoreticalPeakMass -= AminoAcids.getWeightMono(acidSequence[i]);
			out.add((int) Math.round(theoreticalPeakMass));
		}
		return out;
	}
	
	public static ArrayList<Integer> getBIons(Peptide peptide) {
		ArrayList<Integer> out = new ArrayList<Integer>();
		byte [] acidSequence = peptide.getAcidSequence();
		/* b-ion  */
		double theoreticalPeakMass = Properties.leftIonDifference;
		for (int i = 0; i < peptide.getLengthMinusOne(); i++) {
			theoreticalPeakMass += AminoAcids.getWeightMono(acidSequence[i]);
			out.add((int) Math.round(theoreticalPeakMass));
		}
		return out;
	}

}
