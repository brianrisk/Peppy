package Peppy;
import java.io.*;
import java.util.*;

import Utilities.U;


/**
 * The purpose of this class is to have a series of
 * quality control tests.  If any changes are made to any of the
 * scoring mechanisms of JavaGFS, then these tests should
 * be run.
 * @author Brian Risk
 *
 */
public class ReliabilityTester {
	
	public static void runTheGamut(String species) {
		U.startStopwatch();
		ArrayList<Spectrum> spectra = Spectrum.loadSpectraFromFolder("tests/" + species + "/spectra");
		ArrayList<Peptide> peptides = loadHighScoringPeptides(species);
		ArrayList<Peptide> correctPeptides = new  ArrayList<Peptide>();
		ArrayList<String> correctPeptideNames = new ArrayList<String>();
		loadCorrectPeptides(species,correctPeptides,correctPeptideNames );
//		double increment = 0.02;
		double increment = 1.0;
		double thresholdMin = 0.4;
		double thresholdMax = 0.7;
		double powerMin = 0.0;
		double powerMax = 0.5;
		int score;
		int bestScore = 0;
		double bestPower = 0.0;
		double bestThreshold = 0.0;
//		for (double threshold = thresholdMin; threshold < thresholdMax; threshold += increment) {
//			U.p(threshold);
//			for (double power = powerMin; power < powerMax; power += increment) {
//				Properties.peakDifferenceThreshold = threshold;
//				Properties.peakIntensityExponent = power;
//				score = getNumberOfTopRankingMatches("human", spectra, peptides, correctPeptides, correctPeptideNames);
//				if (score > bestScore) {
//					bestScore = score;
//					bestPower = power;
//					bestThreshold = threshold;
//				}
//			}	
//			U.p("SO FAR: Best Score: " + bestScore + " Best Threshold: " + bestThreshold + " Best Power: " + bestPower);
//			U.printTimeRemaining((threshold - thresholdMin) / (thresholdMax - thresholdMin));
//		}
		

		for (double ionDifference = -20.0 ; ionDifference < 40.0; ionDifference += increment) {
			Properties.yIonDifference = ionDifference;
			score = getNumberOfTopRankingMatches("human", spectra, peptides, correctPeptides, correctPeptideNames);
			if (score > bestScore) {
				bestScore = score;
				bestPower = ionDifference;
			}
			U.p("ionDifference: " + ionDifference + "Score: " + score + ".  SO FAR best ion difference: " + bestPower + " with a score of " + bestScore);
		}	
		U.p("best ion difference: " + bestPower);
		
//		U.p("Best Score: " + bestScore);
//		U.p("Best Threshold: " + bestThreshold);
//		U.p("Best Power: " + bestPower);
		U.stopStopwatch();
	}
	
	public static int getNumberOfTopRankingMatches(String species, ArrayList<Spectrum> spectra, ArrayList<Peptide> peptides, ArrayList<Peptide> correctPeptides, ArrayList<String> correctPeptideNames) {
		ArrayList<SpectrumPeptideMatch> matches = JavaGFS.asynchronousDigestion(peptides, spectra, null);
		int out = 0;
		for (int i = 0; i < correctPeptides.size(); i++) {	
			//We've loaded the string.  Now see if we have it as a match to the given spectrum
			ArrayList<SpectrumPeptideMatch> spectrumMatches = getMatchesWithSpectrumName(correctPeptideNames.get(i), matches);
			if (spectrumMatches.size() == 0) continue;
			//We're sorting so that if the peptide is in our list of matches, we know how highly it is ranked.
			Collections.sort(spectrumMatches);
			SpectrumPeptideMatch topMatch = spectrumMatches.get(0);
			if (topMatch.getPeptide().getSequence().equals(correctPeptides.get(i).getSequence())) {
				out++;
			}
		}
		return out;
	}
	
	public static void loadCorrectPeptides(String species, ArrayList<Peptide> correctPeptides, ArrayList<String> correctPeptideNames) {
		//go through each file in our peptides folder
		File peptideFolder = new File("tests/" + species + "/peptides");
		File [] peptideFiles = peptideFolder.listFiles();
		
		for (int peptideFileIndex = 0; peptideFileIndex < peptideFiles.length; peptideFileIndex++) {
			//only want visible files
			if (peptideFiles[peptideFileIndex].isHidden()) continue;
			
			//only want valid file types
			String fileName = peptideFiles[peptideFileIndex].getName();
			if (!fileName.endsWith(".dta") && !fileName.endsWith(".pkl") && !fileName.endsWith(".txt")) continue;
			
			//load in our string line
			String peptideString = "";
			try {
				BufferedReader br = new BufferedReader(new FileReader(peptideFiles[peptideFileIndex]));
				//read the first line;
				peptideString = br.readLine();
				//close;
				br.close();
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
			//Okay, we've got a valid peptide file
			if (peptideString == null) continue;
			peptideString = peptideString.trim();
			if (peptideString.equals("")) continue;
			
			//make sure this peptide has a corresponding spectrum file
			File spectrumFile = new File("tests/" + species + "/spectra/" + peptideFiles[peptideFileIndex].getName());
			if (!spectrumFile.exists()) continue;
			correctPeptides.add(new Peptide(peptideString));
			correctPeptideNames.add(peptideFiles[peptideFileIndex].getName());
			
		}
	}
	
	
	public static void testReliability(String species) {
		U.p("Running reliability test.");
		
		//Load our spectra
		ArrayList<Spectrum> spectra = Spectrum.loadSpectraFromFolder("tests/" + species + "/spectra");
		U.p("loaded " +spectra.size() + " spectra.");
		
		
		ArrayList<Peptide> peptides = loadHighScoringPeptides(species);
		
		//loading peptides from a protein database
//		ArrayList<Peptide> peptides = ProteinDigestion.getPeptidesFromProteinFile(new File("tests/databases/uniprot_sprot.fasta"));
		
		ArrayList<SpectrumPeptideMatch> matches = JavaGFS.asynchronousDigestion(peptides, spectra, null);
		
		//Initialize HMM Score
		U.p("testing HMM here!");
		HMMScore.HMMClass.HmmSetUp();
		for (SpectrumPeptideMatch match: matches) {
			match.calculateHMM();
		}
		SpectrumPeptideMatch.setSortParameter(SpectrumPeptideMatch.SORT_BY_HMM);
		
		Collections.sort(matches);
		for (SpectrumPeptideMatch match: matches) {
			U.p(match.getScoreHMM());
		}
		
		//go through each file in our peptides folder
		File peptideFolder = new File("tests/" + species + "/peptides");
		File [] peptideFiles = peptideFolder.listFiles();
		int numberFound = 0;
		int rankTotal = 0;
		
		int[] rankTotals = new int[Properties.maximumNumberOfMatchesForASpectrum];
		
		
		for (int peptideFileIndex = 0; peptideFileIndex < peptideFiles.length; peptideFileIndex++) {
			//only want visible files
			if (peptideFiles[peptideFileIndex].isHidden()) continue;
			
			//only want valid file types
			String fileName = peptideFiles[peptideFileIndex].getName();
			if (!fileName.endsWith(".dta") && !fileName.endsWith(".pkl") && !fileName.endsWith(".txt")) continue;
			
			//load in our string line
			String peptideString = "";
			try {
				BufferedReader br = new BufferedReader(new FileReader(peptideFiles[peptideFileIndex]));
				//read the first line;
				peptideString = br.readLine();
				//close;
				br.close();
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
			//Okay, we've got a valid peptide file
			if (peptideString == null) continue;
			peptideString = peptideString.trim();
			if (peptideString.equals("")) continue;
			
			//make sure this peptide has a corresponding spectrum file
			File spectrumFile = new File("tests/" + species + "/spectra/" + peptideFiles[peptideFileIndex].getName());
			if (!spectrumFile.exists()) continue;
			
			/*
			 * Here we include the peptide we know to be correct and see how it would score	
			 */
			Peptide peptide = new Peptide(peptideString);
			//This assumes only one spectrum per file, but that should be the case with these test cases.
			Spectrum spectrum = Spectrum.loadSpectra(spectrumFile).get(0);
			Sequence sequence = new Sequence(new File("not important")); 
			SpectrumPeptideMatch realMatch = new SpectrumPeptideMatch(spectrum, peptide, sequence);
//			realMatch.calculateHMM();
			matches.add(realMatch);
			

			//We've loaded the string.  Now see if we have it as a match to the given spectrum
			ArrayList<SpectrumPeptideMatch> spectrumMatches = getMatchesWithSpectrumName(peptideFiles[peptideFileIndex].getName(), matches);
			//We're sorting so that if the peptide is in our list of matches, we know how highly it is ranked.
			Collections.sort(spectrumMatches);
			boolean matchFound = false;
			int maxIndex = Properties.maximumNumberOfMatchesForASpectrum;
			if (spectrumMatches.size() < maxIndex) maxIndex = spectrumMatches.size();
			for (int spectrumIndex = 0; spectrumIndex < maxIndex; spectrumIndex++) {
				SpectrumPeptideMatch topMatch = spectrumMatches.get(spectrumIndex);
				if (topMatch.getPeptide().getSequence().equals(peptideString)) {
					numberFound++;
					rankTotals[spectrumIndex]++;
					matchFound = true;
					rankTotal += spectrumIndex + 1;
					break;
				}
			}
			//If a match isn't found, print out a quick report
//			if (!matchFound) {
//				U.p(
//						spectrumMatches.get(0).getSpectrum().getFile().getName() + "\t" + 
//						spectrumMatches.get(0).getPeptide().getSequence() + " vs " + 
//						peptideString + "\tscore: " + 
//						spectrumMatches.get(0).getScoreMSMSFit() / realMatch.getScoreMSMSFit());
//			}
		}
		U.p("");
		//U.p("peptide count:" + peptideCount);
		U.p("number found in top " + Properties.maximumNumberOfMatchesForASpectrum + ": " + numberFound);
		U.p("precision: " + (double) rankTotals[0] / spectra.size());
		U.p("accuracy: " + (double) numberFound / spectra.size());
		//U.p("average rank: " + (double) rankTotal / numberFound);
		for (int i = 0; i < rankTotals.length; i++) {
			U.p("Rank " + (i + 1) + ": " + rankTotals[i]);
		}
	}
	
	
	private static ArrayList<SpectrumPeptideMatch> getMatchesWithSpectrumName(String spectrumName, ArrayList<SpectrumPeptideMatch> theseMatches) {
		ArrayList<SpectrumPeptideMatch> out = new ArrayList<SpectrumPeptideMatch>();
		for (int i = 0; i < theseMatches.size(); i++) {
			SpectrumPeptideMatch match = theseMatches.get(i);
			if (match.getSpectrum().getFile().getName().equals(spectrumName)) {
				out.add(match);
			}
		}
		return out;
	}
	
	public static void makeSureWeAreProperlyDigestingTheGenome(String species) {
		U.startStopwatch();
		ArrayList<Sequence> sequences = Sequence.loadSequences(new File("tests/" + species + "/sequences"));
		//load peptides to see if we are properly extracting
		ArrayList<Peptide> peptides = sequences.get(0).extractPeptides();
		
		//go through each file in our peptides folder
		File peptideFolder = new File("tests/" + species + "/peptides");
		File [] peptideFiles = peptideFolder.listFiles();
		
		for (int peptideFileIndex = 0; peptideFileIndex < peptideFiles.length; peptideFileIndex++) {
			//only want visible files
			if (peptideFiles[peptideFileIndex].isHidden()) continue;
			
			//only want valid file types
			String fileName = peptideFiles[peptideFileIndex].getName();
			if (!fileName.endsWith(".dta") && !fileName.endsWith(".pkl") && !fileName.endsWith(".txt")) continue;
			
			//load in our string line
			String peptideString = "";
			try {
				BufferedReader br = new BufferedReader(new FileReader(peptideFiles[peptideFileIndex]));
				//read the first line;
				peptideString = br.readLine();
				//close;
				br.close();
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
			//Okay, we've got a valid peptide file
			if (peptideString == null) continue;
			peptideString = peptideString.trim();
			if (peptideString.equals("")) continue;
			
			//make sure this peptide has a corresponding spectrum file
			File spectrumFile = new File("tests/" + species + "/spectra/" + peptideFiles[peptideFileIndex].getName());
			if (!spectrumFile.exists()) continue;
			
			boolean peptideFound = false;
			for (int peptideIndex = 0; peptideIndex < peptides.size(); peptideIndex++) {
				if (peptideString.equals(peptides.get(peptideIndex).getSequence())) {
					peptideFound = true;
					continue;
				}
			}
			//if (!peptideFound) U.p(peptide + " " + peptideFiles[peptideFileIndex].getName());
			if (!peptideFound) U.p(peptideString);
		}
		U.stopStopwatch();
	}
	

	/**
	 * This is a setup method.  What it does is goes through the genome you give it in the
	 * "sequences" directory and finds the highest matching peptides.  It then outputs 
	 * that list of high-scoring peptides to a file.  This file can then be read in
	 * for later scoring comparisons to save time.
	 * @param matches
	 * @param species
	 */
	public static void exportHighScoringPeptides(String species) {
		//Load our spectra
		ArrayList<Spectrum> spectra = Spectrum.loadSpectraFromFolder("tests/" + species + "/spectra");
		U.p("loaded " +spectra.size() + " spectra.");
		
		//Get references to our sequence files -- no nucleotide data is loaded at this point
		ArrayList<Sequence> sequences = Sequence.loadSequences(new File("tests/" + species + "/sequences"));
		
		//initialize our ArrayList of matches
		ArrayList<SpectrumPeptideMatch> matches = new ArrayList<SpectrumPeptideMatch>();

		//loop through each sequence in the sequences ArrayList
		for (int sequenceIndex = 0; sequenceIndex < sequences.size(); sequenceIndex++) {
			Sequence sequence = sequences.get(sequenceIndex);		
			matches.addAll(JavaGFS.asynchronousDigestion(sequence, spectra));
		}
		try {
			//append to existing file
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("tests/" + species + "/highScoringPeptides.txt")));

			//loop through each sequence in the sequences ArrayList
			for (int matchIndex = 0; matchIndex < matches.size(); matchIndex++) {
				SpectrumPeptideMatch match = matches.get(matchIndex);		
				pw.println(match.getPeptide().getSequence());
				double mass1 = match.getPeptide().getMass();
				double mass2 = match.getPeptide().calculateMass();
				if (Math.abs(mass1 - mass2) > 0.0001) U.p(mass2 - mass1);
			}
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void exportPeptidesInCommonWithDatabase(String species) {
		//loading peptides from a protein database
		ArrayList<Peptide> peptides = ProteinDigestion.getPeptidesFromProteinFile(new File("tests/databases/uniprot_sprot.fasta"));
		
		try {
			//Set up our output file
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("tests/" + species + "/peptidesInCommonWithDatabase.txt")));
			
			//go through each file in our peptides folder
			File peptideFolder = new File("tests/" + species + "/peptides");
			File [] peptideFiles = peptideFolder.listFiles();
			
			for (int peptideFileIndex = 0; peptideFileIndex < peptideFiles.length; peptideFileIndex++) {
				//only want visible files
				if (peptideFiles[peptideFileIndex].isHidden()) continue;
				
				//only want valid file types
				String fileName = peptideFiles[peptideFileIndex].getName();
				if (!fileName.endsWith(".dta") && !fileName.endsWith(".pkl") && !fileName.endsWith(".txt")) continue;
				
				//load in our string line
				String peptideString = "";
				try {
					BufferedReader br = new BufferedReader(new FileReader(peptideFiles[peptideFileIndex]));
					//read the first line;
					peptideString = br.readLine();
					//close;
					br.close();
				} catch (FileNotFoundException e) {
					e.printStackTrace();
				} catch (IOException e) {
					e.printStackTrace();
				}
				
				//Okay, we've got a valid peptide file
				if (peptideString == null) continue;
				peptideString = peptideString.trim();
				if (peptideString.equals("")) continue;
				
				//make sure this peptide has a corresponding spectrum file
				File spectrumFile = new File("tests/" + species + "/spectra/" + peptideFiles[peptideFileIndex].getName());
				if (!spectrumFile.exists()) continue;
				
				Peptide peptide = new Peptide(peptideString);
				int peptideIndex = ScoringThread.findFirstIndexWithGreaterMass(peptides, peptide.getMass() - .01);
				double peptideMassButBigger = peptide.getMass() + .01;
				for (int i = peptideIndex; i < peptides.size(); i++) {
					if (peptide.getSequence().equals(peptides.get(i).getSequence())) {
						pw.println(peptide.getSequence());
						break;
					}
					if (peptides.get(i).getMass() > peptideMassButBigger) {
						break;
					}
				}
				
			}
			pw.flush();
			pw.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private static ArrayList<Peptide> loadHighScoringPeptides(String species) {
		ArrayList<Peptide> peptides = new ArrayList<Peptide>();
		try {
			BufferedReader br = new BufferedReader(new FileReader("tests/" + species + "/highScoringPeptides.txt"));
			String line = br.readLine();
			while (line != null) {
				peptides.add(new Peptide(line));
				line = br.readLine();
			}
			br.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		Collections.sort(peptides);
		return peptides;
	}
}
 