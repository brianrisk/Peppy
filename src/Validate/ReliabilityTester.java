package Validate;
import java.io.*;
import java.util.*;

import Peppy.Peppy;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.ProteinDigestion;
import Peppy.ScoringThreadServer;
import Peppy.ScoringThread;
import Peppy.Sequence;
import Peppy.Spectrum;
import Peppy.Match;
import Utilities.U;


/**
 * CAUTION: not designed to be user friendly
 * 
 * The purpose of this class is to have a series of
 * quality control tests.  If any changes are made to any of the
 * scoring mechanisms of JavaGFS, then these tests should
 * be run.
 * @author Brian Risk
 *
 */
public class ReliabilityTester {
	
	public static void main(String slwen[]) {
		U.p("testing reliablity.");
//		ReliabilityTester.runTheGamut("USP");
		ReliabilityTester.runTheGamut2("human");
//		ReliabilityTester.testReliability("aurum");
//		ReliabilityTester.testReliability("human");
//		ReliabilityTester.testReliability("ecoli");
//		ReliabilityTester.testReliability("USP");
//		ReliabilityTester.exportHighScoringPeptidesFromSwissProt("USP");
//		ReliabilityTester.exportPeptidesInCommonWithDatabase("human");
//		ReliabilityTester.exportPeptidesInCommonWithDatabase("ecoli");
//		ReliabilityTester.makeSureWeAreProperlyDigestingTheGenome("ecoli");
		U.p("done!");
	}
	
	public static void runTheGamut2(String species) {
		U.p("running the gamut!");
		U.startStopwatch();
		ArrayList<Spectrum> spectra = Spectrum.loadSpectraFromFolder("/Users/risk2/PeppyOverflow/tests/" + species + "/spectra");
		ArrayList<Peptide> peptides = loadHighScoringPeptides(species);
		ArrayList<Peptide> correctPeptides = new  ArrayList<Peptide>();
		ArrayList<String> correctPeptideNames = new ArrayList<String>();
		loadCorrectPeptides(species,correctPeptides,correctPeptideNames );
		double increment = 0.1;
		double thresholdMin = 0.7;
		double thresholdMax = 1.72;
		
//		best YBtrue:1.0999999999999999
//		best YBfalse:1.2
//		best BYtrue:1.3
//		best BYtrue:0.8999999999999999
		
		double YBtrue, YBfalse, BYtrue, BYfalse;
		double bestYBtrue = 0, bestYBfalse = 0, bestBYtrue = 0, bestBYfalse = 0;
		int score;
		int bestScore = 0;
		for (YBtrue = 1.1; YBtrue <= thresholdMax; YBtrue += increment) {
			for (YBfalse = 1.2; YBfalse <= thresholdMax; YBfalse += increment) {
				for (BYtrue = 0.7; BYtrue <= thresholdMax; BYtrue += increment) {
					for (BYfalse = 0.6; BYfalse <= 1; BYfalse += increment) {
						Properties.YBtrue = YBtrue;
						Properties.YBfalse = YBfalse;
						Properties.BYtrue = BYtrue;
						Properties.BYfalse = BYfalse;
						score = getNumberOfTopRankingMatches(species, spectra, peptides, correctPeptides, correctPeptideNames);
						U.p ("score: " + score + "; " + YBtrue + ", " + YBfalse + ", " + BYtrue + ", " + BYfalse);
						if (score > bestScore) {
							bestScore = score;
							bestYBtrue = YBtrue;
							bestYBfalse = YBfalse;
							bestBYtrue = BYtrue;
							bestBYfalse = BYfalse;
							
							U.p();
							U.p("New best scores:");
							U.p("best YBtrue:" + bestYBtrue);
							U.p("best YBfalse:" + bestYBfalse);
							U.p("best BYtrue:" + bestBYtrue);
							U.p("best BYtrue:" + bestBYfalse);
						}
					}
				}
			}
			
		}

		U.p();
		U.p("BEST SCORES:");
		U.p("best YBtrue:" + bestYBtrue);
		U.p("best YBfalse:" + bestYBfalse);
		U.p("best BYtrue:" + bestBYtrue);
		U.p("best BYtrue:" + bestBYfalse);
		U.stopStopwatch();
	}
	
	public static void runTheGamut(String species) {
		U.startStopwatch();
		ArrayList<Spectrum> spectra = Spectrum.loadSpectraFromFolder("/Users/risk2/PeppyOverflow/tests/" + species + "/spectra");
		ArrayList<Peptide> peptides = loadHighScoringPeptides(species);
		ArrayList<Peptide> correctPeptides = new  ArrayList<Peptide>();
		ArrayList<String> correctPeptideNames = new ArrayList<String>();
		loadCorrectPeptides(species,correctPeptides,correctPeptideNames );
		double increment = 0.05;
		double thresholdMin = 0.2;
		double thresholdMax = 0.7;
		double powerMin = 0.2;
		double powerMax = 0.5;
		int score;
		int bestScore = 0;
		double bestPower = 0.0;
		double bestThreshold = 0.0;
		boolean bestHighIntensityCleaning = false;
		boolean bestLocalMaximaCleaning = false;
		for (double threshold = thresholdMin; threshold < thresholdMax; threshold += increment) {
			U.p(threshold);
			for (double power = powerMin; power < powerMax; power += increment) {
				for (int boolcount = 0; boolcount < 2; boolcount++) {
					for (int boolcount2 = 0; boolcount2 < 2; boolcount2++) {
						Properties.highIntensityCleaning = (boolcount == 0);
						Properties.localMaximaCleaning = (boolcount2 == 0);
						Properties.peakDifferenceThreshold = threshold;
						Properties.peakIntensityExponent = power;
						score = getNumberOfTopRankingMatches(species, spectra, peptides, correctPeptides, correctPeptideNames);
						if (score > bestScore) {
							bestScore = score;
							bestPower = power;
							bestThreshold = threshold;
							bestHighIntensityCleaning = (boolcount == 0);
							bestLocalMaximaCleaning = (boolcount2 == 0);
						}
					}
				}
			}	
			U.p("SO FAR: Best Score: " + bestScore + " Best Threshold: " + bestThreshold + " Best Power: " + bestPower);
			U.p("Cleaning: " + bestHighIntensityCleaning + ", " + bestLocalMaximaCleaning);
			U.p("Percent remaining is: " + (threshold - thresholdMin)  * 100 / (thresholdMax - thresholdMin));
		}
		U.p("Best Score: " + bestScore);
		U.p("Best Threshold: " + bestThreshold);
		U.p("Best Power: " + bestPower);
		U.p("Best bestHighIntensityCleaning:" + bestHighIntensityCleaning);
		U.p("Best bestLocalMaximaCleaning:" + bestLocalMaximaCleaning);
		U.stopStopwatch();
	}
	
	public static int getNumberOfTopRankingMatches(String species, ArrayList<Spectrum> spectra, ArrayList<Peptide> peptides, ArrayList<Peptide> correctPeptides, ArrayList<String> correctPeptideNames) {
		ArrayList<Match> matches = (new ScoringThreadServer(peptides, spectra, null)).getMatches();
		int out = 0;
		for (int i = 0; i < correctPeptides.size(); i++) {	
			//We've loaded the string.  Now see if we have it as a match to the given spectrum
			ArrayList<Match> spectrumMatches = getMatchesWithSpectrumName(correctPeptideNames.get(i), matches);
			if (spectrumMatches.size() == 0) continue;
			//We're sorting so that if the peptide is in our list of matches, we know how highly it is ranked.
			Collections.sort(spectrumMatches);
			Match topMatch = spectrumMatches.get(0);
			if (topMatch.getPeptide().getAcidSequence().equals(correctPeptides.get(i).getAcidSequence())) {
				out++;
			}
		}
		return out;
	}
	
	public static void loadCorrectPeptides(String species, ArrayList<Peptide> correctPeptides, ArrayList<String> correctPeptideNames) {
		//go through each file in our peptides folder
		File peptideFolder = new File("/Users/risk2/PeppyOverflow/tests/" + species + "/peptides");
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
			File spectrumFile = new File("/Users/risk2/PeppyOverflow/tests/" + species + "/spectra/" + peptideFiles[peptideFileIndex].getName());
			if (!spectrumFile.exists()) continue;
			correctPeptides.add(new Peptide(peptideString));
			correctPeptideNames.add(peptideFiles[peptideFileIndex].getName());
			
		}
	}
	
	
	public static void testReliability(String species) {
		U.p("Running reliability test.");
		
		//Load our spectra
		ArrayList<Spectrum> spectra = Spectrum.loadSpectraFromFolder("/Users/risk2/PeppyOverflow/tests/" + species + "/spectra");
		U.p("loaded " +spectra.size() + " spectra.");
		
		
		ArrayList<Peptide> peptides = loadHighScoringPeptides(species);
		
		//loading peptides from a protein database
//		ArrayList<Peptide> peptides = ProteinDigestion.getPeptidesFromProteinFile(new File("tests/databases/uniprot_sprot.fasta"));
		
		ArrayList<Match> matches = (new ScoringThreadServer(peptides, spectra, null)).getMatches();
		
		//Initialize HMM Score
//		U.p("testing HMM here!");
//		HMMScore.HMMClass.HmmSetUp();
//		for (SpectrumPeptideMatch match: matches) {
//			match.calculateHMM();
//		}
//		SpectrumPeptideMatch.setSortParameter(SpectrumPeptideMatch.SORT_BY_HMM);
//		
//		Collections.sort(matches);
//		for (SpectrumPeptideMatch match: matches) {
//			U.p(match.getScoreHMM());
//		}
		
		//go through each file in our peptides folder
		File peptideFolder = new File("/Users/risk2/PeppyOverflow/tests/" + species + "/peptides");
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
			File spectrumFile = new File("/Users/risk2/PeppyOverflow/tests/" + species + "/spectra/" + peptideFiles[peptideFileIndex].getName());
			if (!spectrumFile.exists()) continue;
			
			/*
			 * Here we include the peptide we know to be correct and see how it would score	
			 */
			Peptide peptide = new Peptide(peptideString);
			//This assumes only one spectrum per file, but that should be the case with these test cases.
			Spectrum spectrum = Spectrum.loadSpectra(spectrumFile).get(0);
			Sequence sequence = new Sequence(new File("not important")); 
			Match realMatch = new Match(spectrum, peptide, sequence);
//			realMatch.calculateHMM();
			matches.add(realMatch);
			

			//We've loaded the string.  Now see if we have it as a match to the given spectrum
			ArrayList<Match> spectrumMatches = getMatchesWithSpectrumName(peptideFiles[peptideFileIndex].getName(), matches);
			//We're sorting so that if the peptide is in our list of matches, we know how highly it is ranked.
			Collections.sort(spectrumMatches);
			boolean matchFound = false;
			int maxIndex = Properties.maximumNumberOfMatchesForASpectrum;
			if (spectrumMatches.size() < maxIndex) maxIndex = spectrumMatches.size();
			for (int spectrumIndex = 0; spectrumIndex < maxIndex; spectrumIndex++) {
				Match topMatch = spectrumMatches.get(spectrumIndex);
				if (topMatch.getPeptide().getAcidSequence().equals(peptideString)) {
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
	
	
	private static ArrayList<Match> getMatchesWithSpectrumName(String spectrumName, ArrayList<Match> theseMatches) {
		ArrayList<Match> out = new ArrayList<Match>();
		for (int i = 0; i < theseMatches.size(); i++) {
			Match match = theseMatches.get(i);
			if (match.getSpectrum().getFile().getName().equals(spectrumName)) {
				out.add(match);
			}
		}
		return out;
	}
	
	public static void makeSureWeAreProperlyDigestingTheGenome(String species) {
		U.startStopwatch();
		ArrayList<Sequence> sequences = Sequence.loadSequences(new File("/Users/risk2/PeppyOverflow/tests/" + species + "/sequences"));
		//load peptides to see if we are properly extracting
		ArrayList<Peptide> peptides = sequences.get(0).extractMorePeptides();
		
		//go through each file in our peptides folder
		File peptideFolder = new File("/Users/risk2/PeppyOverflow/tests/" + species + "/peptides");
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
			File spectrumFile = new File("/Users/risk2/PeppyOverflow/tests/" + species + "/spectra/" + peptideFiles[peptideFileIndex].getName());
			if (!spectrumFile.exists()) continue;
			
			boolean peptideFound = false;
			for (int peptideIndex = 0; peptideIndex < peptides.size(); peptideIndex++) {
				if (peptideString.equals(peptides.get(peptideIndex).getAcidSequence())) {
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
		U.p("finding high scoring peptides");
		//Load our spectra
		ArrayList<Spectrum> spectra = Spectrum.loadSpectraFromFolder("/Users/risk2/PeppyOverflow/tests/" + species + "/spectra");
		U.p("loaded " +spectra.size() + " spectra.");
		
		//Get references to our sequence files -- no nucleotide data is loaded at this point
		ArrayList<Sequence> sequences = Sequence.loadSequences(new File("/Users/risk2/PeppyOverflow/tests/" + species + "/sequences"));
		
		//initialize our ArrayList of matches
		ArrayList<Match> matches = new ArrayList<Match>();

		//loop through each sequence in the sequences ArrayList
		for (int sequenceIndex = 0; sequenceIndex < sequences.size(); sequenceIndex++) {
			Sequence sequence = sequences.get(sequenceIndex);		
			matches.addAll(Peppy.getMatches(sequence, spectra));
		}
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("/Users/risk2/PeppyOverflow/tests/" + species + "/highScoringPeptides.txt")));
			//loop through each sequence in the sequences ArrayList
			for (int matchIndex = 0; matchIndex < matches.size(); matchIndex++) {
				Match match = matches.get(matchIndex);		
				pw.println(match.getPeptide().getAcidSequence());
			}
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void exportHighScoringPeptidesFromSwissProt(String species) {
		U.p("finding high scoring peptides from swiss prot");
		//Load our spectra
		ArrayList<Spectrum> spectra = Spectrum.loadSpectraFromFolder("/Users/risk2/PeppyOverflow/tests/" + species + "/spectra");
		U.p("loaded " +spectra.size() + " spectra.");
		
		//load SwissProt
		ArrayList<Peptide> peptides = ProteinDigestion.getPeptidesFromFASTA(new File("/Users/risk2/PeppyOverflow/tests/databases/uniprot_sprot.fasta"));
		
		//keep the top 20 matches for each spectrum
		Properties.maximumNumberOfMatchesForASpectrum = 20;
		ArrayList<Match> matches = (new ScoringThreadServer(peptides, spectra, null)).getMatches();
		
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("/Users/risk2/PeppyOverflow/tests/" + species + "/highScoringPeptides.txt")));

			//loop through each sequence in the sequences ArrayList
			for (Match match: matches) {	
				pw.println(match.getPeptide().getAcidSequence());
			}
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void exportPeptidesInCommonWithDatabase(String species) {
		//loading peptides from a protein database
		ArrayList<Peptide> peptides = ProteinDigestion.getPeptidesFromFASTA(new File("/Users/risk2/PeppyOverflow/tests/databases/uniprot_sprot.fasta"));
		
		try {
			//Set up our output file
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("/Users/risk2/PeppyOverflow/tests/" + species + "/peptidesInCommonWithDatabase.txt")));
			
			//go through each file in our peptides folder
			File peptideFolder = new File("/Users/risk2/PeppyOverflow/tests/" + species + "/peptides");
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
				File spectrumFile = new File("/Users/risk2/PeppyOverflow/tests/" + species + "/spectra/" + peptideFiles[peptideFileIndex].getName());
				if (!spectrumFile.exists()) continue;
				
				Peptide peptide = new Peptide(peptideString);
				int peptideIndex = ScoringThread.findFirstIndexWithGreaterMass(peptides, peptide.getMass() - .01);
				double peptideMassButBigger = peptide.getMass() + .01;
				for (int i = peptideIndex; i < peptides.size(); i++) {
					if (peptide.getAcidSequence().equals(peptides.get(i).getAcidSequence())) {
						pw.println(peptide.getAcidSequence());
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
	
	public static ArrayList<Peptide> loadHighScoringPeptides(String species) {
		File file = new File("/Users/risk2/PeppyOverflow/tests/" + species + "/highScoringPeptides.txt");
		return loadHighScoringPeptides(file);
	}
	
	public static ArrayList<Peptide> loadHighScoringPeptides(File file) {
		ArrayList<Peptide> peptides = new ArrayList<Peptide>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
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
 