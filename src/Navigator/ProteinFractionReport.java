package Navigator;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

import Peppy.Protein;
import Peppy.ProteinCoverage;
import Peppy.Sequence_Protein;
import Peppy.U;

/**
 * Reports specific for techniques used in Hanash lab at MD Anderson
 * 
 * 
 * @author Brian Risk
 *
 */
public class ProteinFractionReport {
	

	public static void main(String [] args) {
		/* load protein database */
		Hashtable<String, ProteinCoverage> proteins = new Hashtable<String,ProteinCoverage>();
		File proteinFile = new File("/Users/risk2/PeppyData/public/sequences/protein/UniProt_Human_2012_03.fasta");
		Sequence_Protein proteinDatabase = new Sequence_Protein(proteinFile);
		ArrayList<Protein> proteinsLoaded = proteinDatabase.getProteinsFromDatabase(false, false);	
		for (Protein protein: proteinsLoaded) {	
			proteins.put(protein.getName(), new ProteinCoverage(protein));
		}
		
		/* load matches */
		ArrayList<Match> matches = new ArrayList<Match>();
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/Documents/MD Anderson/Visit results/group1/1 group1 - UniProt_Human_2012_03.fasta/report.txt")));
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/Documents/MD Anderson/Visit results/group2/1 group2 - UniProt_Human_2012_03.fasta/report.txt")));
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/Documents/MD Anderson/Visit results/group3/1 group3 - UniProt_Human_2012_03.fasta/report.txt")));
		matches.addAll(Match.loadMatches(new File("/Users/risk2/Documents/MD Anderson/Visit results/group4/1 group4 - UniProt_Human_2012_03.fasta/report.txt")));
		matches.addAll(Match.loadMatches(new File("/Users/risk2/Documents/MD Anderson/Visit results/group5/1 group5 - UniProt_Human_2012_03.fasta/report.txt")));
		matches.addAll(Match.loadMatches(new File("/Users/risk2/Documents/MD Anderson/Visit results/group6/1 group6 - UniProt_Human_2012_03.fasta/report.txt")));
		
		/* see which peptides are found in only one protein form */
		Hashtable<String, Boolean> spectraWithOnePeptideMatch = new Hashtable<String, Boolean>();
		String spectrumMD5;
		for (Match match: matches) {
			spectrumMD5 = match.getString("spectrumMD5");
			Boolean found = spectraWithOnePeptideMatch.get(spectrumMD5);
			if (found == null) {
				spectraWithOnePeptideMatch.put(spectrumMD5, true);
			} else {
				spectraWithOnePeptideMatch.put(spectrumMD5, false);
			}
		}
		
		/* all the fraction names */
		ArrayList<String> fractionLabels = new ArrayList<String>();
		fractionLabels.add("01to21");
		fractionLabels.add("22to24");
		fractionLabels.add("25to26");
		fractionLabels.add("27to28");
		fractionLabels.add("29to30");
		fractionLabels.add("31to32");
		fractionLabels.add("33to34");
		fractionLabels.add("35to36");
		fractionLabels.add("37to38");
		fractionLabels.add("39to41");
		fractionLabels.add("42to43");
		fractionLabels.add("44to45");
		fractionLabels.add("46to47");
		fractionLabels.add("48to49");
		fractionLabels.add("50to51");
		fractionLabels.add("52to53");
		fractionLabels.add("54to55");
		fractionLabels.add("56to57");
		fractionLabels.add("58to59");
		fractionLabels.add("60to62");
		fractionLabels.add("63to65");
		fractionLabels.add("66to68");
		fractionLabels.add("69to77");
		fractionLabels.add("78to84");
		
		
		/* associating fractions with integer values so we can see which are far apart */
		Hashtable<String, Integer> fractionValues = new Hashtable<String, Integer>();
		for (int i = 0; i < fractionLabels.size(); i++) {
			fractionValues.put(fractionLabels.get(i), i + 1);
		}
		
		
		Hashtable<String, ArrayList<Integer>> proteinsFoundWithUniquePeptides = new Hashtable<String, ArrayList<Integer>>();
		Hashtable<String, String> uniquePeptideProteinNames = new Hashtable<String, String>();
		for (Match match: matches) {
			spectrumMD5 = match.getString("spectrumMD5");
			if (spectraWithOnePeptideMatch.get(spectrumMD5)) {
				String peptideSequence = match.getString("peptideSequence");;
				String fractionName = getFractionFromString(match.getFile("FilePath").getName());
				String proteinName = match.getString("SequenceName");
				
				/* adding to our list of which identifying peptides map to which proteins */
				uniquePeptideProteinNames.put(peptideSequence, proteinName);
				
				
				
				ArrayList<Integer> foundFractionIntegers = proteinsFoundWithUniquePeptides.get(proteinName);
				if (foundFractionIntegers == null) foundFractionIntegers = new ArrayList<Integer>();
				foundFractionIntegers.add(fractionValues.get(fractionName));
				proteinsFoundWithUniquePeptides.put(proteinName, foundFractionIntegers);
				
			}
			
		}
		
		ArrayList<String> listOfFoundProteins = new ArrayList<String>(proteinsFoundWithUniquePeptides.keySet());
		for (String proteinName: listOfFoundProteins) {
			ArrayList<Integer> foundFractionIntegers = proteinsFoundWithUniquePeptides.get(proteinName);
			Collections.sort(foundFractionIntegers);
			int previous = -1;
			boolean possibleInterest = false;
			for (int fractionInteger: foundFractionIntegers) {
				if (previous == -1) previous = fractionInteger;
				if (fractionInteger - previous > 7) {
					possibleInterest = true;
				}
				previous = fractionInteger;
			}
			
			if (possibleInterest) {
				U.p("proteinsOfInterest.put(\"" + proteinName + "\", \"" + proteinName + "\");");
				
				Hashtable<Integer, Integer> reducedFractions = new Hashtable<Integer, Integer>();
				for (int fractionInteger: foundFractionIntegers) {
					reducedFractions.put(fractionInteger, fractionInteger);
				}
				
				StringBuffer println = new StringBuffer();
				ArrayList<Integer> rfValues = new ArrayList<Integer>(reducedFractions.values());
				Collections.sort(rfValues);
				for (int fractionInteger: rfValues) {
					println.append(fractionInteger + ", ");
				}
//				U.p(println);
//				U.p();
				
			}
		}
			
		
	}
	
	
	public static String getFractionFromString(String string) {
		return string.substring(string.indexOf("_SG") + 3, string.length() - 4);
	}

}
