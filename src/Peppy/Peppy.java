package Peppy;
import java.io.File;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;
import java.util.ArrayList;
import java.util.Collections;

import Reports.HTMLReporter;
import Reports.TextReporter;
import Utilities.U;

/**
 * Peppy
 * Designed with the following goals:
 * 1) More simple code to promote open source development
 * 2) Easy proteogenomic mapping
 * 3) better multi-threading
 * @author Brian Risk
 *
 */
public class Peppy {
	
	/* So that we may report the total amount of peptides found */
	static int peptideTally = 0;
	
	/* track how much memory we have used */
	static MemoryUsage memoryUsage;
	static long maxMemoryUsed = 0;
	
	
	public static void main(String [] args) {
		/* set up initial state */
		init(args);
		
		/*  hello! */
		printGreeting();
		
		/* do the work */
		runJobs(args);
		
		/* i'm finished! */
		printFarewell();
	}
	
	
	public static void runPeppy(String [] args) {
		U.startStopwatch();
		peptideTally = 0;
		
		//create new report directory
		File reportDir = new File(Properties.reportDirectory, Properties.spectraDirectoryOrFile.getName() + "_" + System.currentTimeMillis());
		
		//save our properties
		Properties.generatePropertiesFile(reportDir);
		
		//Load our spectra
		U.p("loading spectra...");
		ArrayList<Spectrum> spectra = Spectrum.loadSpectra();
		U.p("loaded " +spectra.size() + " spectra.");
		
		//Get references to our sequence files -- no nucleotide data is loaded at this point
		ArrayList<Sequence> sequences = Sequence.loadSequenceFiles(Properties.sequenceDirectoryOrFile);
		
		//initialize our ArrayList of matches
		ArrayList<Match> matches = null;
		
		if (Properties.useSpliceVariants) {
			//gets the first nucleotide sequence in the first sequence file
			Sequence_DNA sequenceFile = (Sequence_DNA) sequences.get(0);
			RNA_Sequence rna = new RNA_Sequence(sequenceFile, sequenceFile.getNucleotideSequences().get(0), Properties.sequenceRegionStart, Properties.sequenceRegionStop);
			
			U.p("digesting with splice variants...");
			RNA_Digestor rnaDigestor = new RNA_Digestor(rna);
			ArrayList<Peptide> peptides  = rnaDigestor.getPeptides();
			U.p("generated " + peptides.size() + " peptides");
			
			U.p("getting matches...");
			//TODO get rid of this getMatches function when this is overhauled
			matches = getMatchesWithPeptides(peptides, spectra);
			
		} else {	
			if (Properties.useSequenceRegion) {
				U.p("digesting part of sequence");
				ArrayList<Sequence> oneSequenceList = new ArrayList<Sequence>();
				oneSequenceList.add(sequences.get(0));
				sequences = oneSequenceList;
			}
			
			/* nothing weird.  Just do a normal search */
			matches = getMatches(sequences, spectra);

		}
		
		/*print peptide tally */
		U.p("The total number of peptides produced is: " + peptideTally);
		
		/* regions */
		Regions regions = new Regions(matches, sequences, spectra);
		U.p("we found this many regions: " + regions.getRegions().size());
	
		/* localized PTM search */
		U.p("performing localized PTM search");
		Properties.matchConstructor = new MatchConstructor("Peppy.Match_IMP_VariMod");
		Properties.searchModifications = true;
		Properties.modificationLowerBound = -100;
		Properties.modificationUpperBound = 100;
		ArrayList<Peptide> peptidesInRegions = getPeptidesInRegions(regions.getRegions(), sequences, false);
		U.p("there are this many peptides in the regions: " + peptidesInRegions.size());
		matches = getMatchesWithPeptides(peptidesInRegions, spectra);
		
		/* rerun the regions analysis */
		regions = new Regions(matches, sequences, spectra);
		
		U.p("creating text reports");
		TextReporter textReport = new TextReporter(matches, spectra, sequences, reportDir);
		textReport.generateFullReport();
		
		if (Properties.createHTMLReport) {
			U.p("creating HTML reports");
			HTMLReporter report = new HTMLReporter(matches, spectra, sequences, reportDir);
			report.generateFullReport();
		}	
		
		/* creating regions report */
		regions.createReport(reportDir);
		
		
		//clear out memory
		regions.clearRegions();
		matches.clear();
		spectra.clear();
		sequences.clear();
		System.gc();
		
		U.stopStopwatch();
	}
	
	public static void runJobs(String [] args) {
		File jobsDir = new File("jobs");
		File[] potentialJobsFiles = jobsDir.listFiles();
		ArrayList<File> jobFiles = new ArrayList<File>();
		if (potentialJobsFiles != null) {
			for (int i = 0; i < potentialJobsFiles.length; i++) {
				if (potentialJobsFiles[i].getName().toLowerCase().endsWith(".txt")) {
					jobFiles.add(potentialJobsFiles[i]);
				}	
			}
		}
		if (jobFiles.size() == 0) {
			U.p("no jobs in jobs folder.  running according to main properties file");
			runPeppy(null);
		} else {
			U.p("running " + jobFiles.size() + " jobs");
			for (int i = 0; i < jobFiles.size(); i++) {
				U.p("running job " + (i + 1) + "; " + jobFiles.get(i).getName());
				init(args);
				init(jobFiles.get(i).getAbsolutePath());
				runPeppy(null);
			}
		}
	}


	public static void init(String propertiesFile) {
		System.setProperty("java.awt.headless", "true"); 
		Properties.loadProperties(propertiesFile);
		AminoAcids.init();
		
		/* track memory */
		MemoryMXBean mbean = ManagementFactory.getMemoryMXBean();
		memoryUsage = mbean.getHeapMemoryUsage();
	}
	
	
	public static void init(String [] args) {
		if (args.length == 0) {
			init("properties.txt");
		} else {
			init(args[0]);
			U.p("Initializing with properties file: " + args[0]);
		}
	}
	
	
	public static ArrayList<Match> getReverseMatches(ArrayList<Sequence> sequences, ArrayList<Spectrum> spectra) {
		return getMatches(sequences, spectra, true);
	}
	
	
	/**
	 * Assumes that we are doing the normal, forwards digestion of our sequences
	 * @param sequences
	 * @param spectra
	 * @return
	 */
	public static ArrayList<Match> getMatches(ArrayList<Sequence> sequences, ArrayList<Spectrum> spectra) {
		if (Properties.useIsotopeLabeling) {
			return getLabeledMatches(sequences, spectra);
		} else {
		return getMatches(sequences, spectra, Properties.useReverseDatabase);
		}
	}

	/**
	 * This is the heart of Peppy where the grand symphony takes place.
	 * 
	 * @param sequences our list of sequences where we will be getting our peptides
	 * @param spectra our list of spectra
	 * @param isReverse if we are doing a normal, forwards search or if this is a null, reverse search
	 * @return
	 */
	private static ArrayList<Match> getMatches(ArrayList<Sequence> sequences, ArrayList<Spectrum> spectra, boolean isReverse) {
		
		/* define the chunk of spectra on which we will work */
		int spectraStart = 0;
		int spectraStop = Properties.numberOfSpectraPerSegment;
		if (spectraStop > spectra.size()) spectraStop = spectra.size();
		
		/* initialize matches and set the initial capacity to be a multiple of the spectra size */
		ArrayList<Match> matches = new ArrayList<Match>(spectra.size() * Properties.maximumNumberOfMatchesForASpectrum);
		
		/* This loop breaks when we have run our last chunk of spectra */
		while (true) {
			U.p("working on spectra " + spectraStart + " to " + spectraStop);
			ArrayList<Match> segmentMatches = new ArrayList<Match>();
			
			/* get a manageable segment of spectra */
			ArrayList<Spectrum> spectraSegment = new ArrayList<Spectrum>();
			for (int i = spectraStart; i < spectraStop; i++) {
				spectraSegment.add(spectra.get(i));
			}
			
			/* track which sequence we are examining */
			int sequenceIndex = 0;
			
				
			/* loops until we have gone through all of our sequences */
			while (true) {
				/* check our memory situation */
				if (maxMemoryUsed < memoryUsage.getUsed()) {
					maxMemoryUsed = memoryUsage.getUsed();
				}
				
				/* Extract a decent size of peptides.  Sequences may be short, so this
				 * goes through each sequence and extracts peptides until a desired
				 * threshold has been been reached.
				 */
				
				/*  Initialize our list of peptides */
				ArrayList<Peptide> peptides = new ArrayList<Peptide>();
				
				/* This is where we get a chunk of peptides */
				ArrayList<Peptide> peptideSegment = new ArrayList<Peptide>();
				
				/* if we use a region, extract that region, else go through all sequences getting a chunk at a time. */
				if (Properties.useSequenceRegion) {
					
					/* collect peptides */
					peptideSegment = sequences.get(sequenceIndex).extractMorePeptides(isReverse);
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
						peptideSegment = sequences.get(sequenceIndex).extractMorePeptides(isReverse);
						
							
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
				
				/* sort our peptide list by mass */
				Collections.sort(peptides);
				
				/* report */
				U.p("we are processing a chunk of peptides this size: " + peptides.size());
				
				/* Where we store our matches for this batch of peptides */
				ArrayList<Match> newMatches;
				
					
				/* keep a count of all peptides */
				peptideTally += peptides.size();
				
				/* run the scoring threads to get the matches from our set of peptides and spectra */
				newMatches = (new ScoringThreadServer(peptides, spectraSegment)).getMatches();
				
				/* add the new matches with tolerable IMP values to the segment matches */
				for (Match match: newMatches) {
					if (match.getIMP() <= Properties.maxIMP) segmentMatches.add(match);
				}
				
				/* free up memory */
				newMatches.clear();
				peptides.clear();
				System.gc();
					
				/* break if we have covered our last sequence or if we are only using a region */
				if (sequenceIndex == sequences.size() || Properties.useSequenceRegion) {
					break;
				}

			}

			
			/* Here we do some basic processing and cleaning of the matches */
			removeDuplicateMatches(segmentMatches);
			assignConfidenceValuesToMatches(segmentMatches, spectra);
			assignRankToMatches(segmentMatches);
			segmentMatches = removePoorMatches(segmentMatches);
			
			
			/* add segment matches to the full list of matches */
			matches.addAll(segmentMatches);
			
			/* clear out memory */
			segmentMatches.clear();
			System.gc();
			
			/* increment our spectrum segment markers */
			if (spectraStop == spectra.size()) break;
			spectraStart = spectraStop;
			spectraStop += Properties.numberOfSpectraPerSegment;
			if (spectraStop > spectra.size()) spectraStop = spectra.size();
			
		}
		
		assignRepeatedPeptideCount(matches);
				
		return matches;
	}
	
	/**
	 * Gets matches where a list of peptides is already derived
	 * @param peptides
	 * @param spectra
	 * @param sequence_DNA
	 * @return
	 */
	public static ArrayList<Match> getMatchesWithPeptides(ArrayList<Peptide> peptides, ArrayList<Spectrum> spectra) {
		peptideTally += peptides.size();
		
		//This is where the bulk of the processing in long jobs takes
		ArrayList<Match> matches  = (new ScoringThreadServer(peptides, spectra)).getMatches();

		//Add only matches with a decent e value	
//		U.p("removing duplicate matches...");
//		removeDuplicateMatches(matches);
		U.p("assigning rank to matches...");
		assignRankToMatches(matches);
		U.p("assigning repeated peptide count to matches...");
		assignRepeatedPeptideCount(matches);
		U.p("assigning confidence values to matches...");
		assignConfidenceValuesToMatches(matches, spectra);
		U.p("removing poor matches...");
		matches = removePoorMatches(matches);
		
		return matches;
	}
	
	public static ArrayList<Match> getLabeledMatches(ArrayList<Sequence> sequences, ArrayList<Spectrum> spectra) {
		U.p("performing isotope labeling analysis");
		
		/* getting labeled matches */
		ArrayList<Match> labeledMatches = getMatches(sequences, spectra, Properties.useReverseDatabase);
		
		/*getting unlabled... */
		Properties.useIsotopeLabeling = false;
		ArrayList<Match> unLabeledMatches = getMatches(sequences, spectra, Properties.useReverseDatabase);
		
		/* resetting the property */
		//TODO we should set up a different property
		Properties.useIsotopeLabeling = true;
		
		/*initializing our final output */
		ArrayList<Match> matches = new ArrayList<Match>(unLabeledMatches.size() + labeledMatches.size());
		
		
		boolean next = false;
		for (Match match: labeledMatches) {
			if (next) {
				next = false;
				continue;
			}
			for (Match other: unLabeledMatches) {
				if (match.getPeptide().equals(other.getPeptide())) {
					next = true;
					match.setScore(match.getScore() + 6);
					match.calculateEValue();
					match.setHasIsotopeConfirmation(true);
					break;
				}
			}			
		}
		
		/* now reverse the order -- do unLabeledMatches first */
		for (Match match: unLabeledMatches) {
			if (next) {
				next = false;
				continue;
			}
			for (Match other: labeledMatches) {
				if (match.getPeptide().equals(other.getPeptide())) {
					next = true;
					match.setScore(match.getScore() + 10);
					match.calculateEValue();
					match.setHasIsotopeConfirmation(true);
					break;
				}
			}	
		}
		
		/* combine sets */
		matches.addAll(labeledMatches);
		matches.addAll(unLabeledMatches);
		
		return matches;
	}


	/**
	 * First, this sets the ID of a match.
	 * 
	 * For every spectrum there is one or a set of #1 ranking matches, #2 ranking matches, etc
	 * There can be a set because the same peptide can appear in multiple places within a genome
	 * @param matches
	 */
	public static void assignRankToMatches(ArrayList<Match> matches) {
		/* first make sure we have matches */
		if (matches.size() == 0) return;
		
		/* sort by spectrum ID, then score */
		Match.setSortParameter(Match.SORT_BY_SPECTRUM_ID_THEN_SCORE);
		Collections.sort(matches);
		Match match = matches.get(0);
		Match previousMatch = matches.get(0);
		//set for the first
		match.rank = 1;
		int rank = 1;
		for (int i = 1; i < matches.size(); i++) {
			//see if these are matches for a different spectrum
			match = matches.get(i);
			
			/* assign an ID to each match */
			match.setId(i);
			
			if (match.getSpectrum().getId() != previousMatch.getSpectrum().getId()) {
				rank = 1;
			} else {
				if (match.getScore() == previousMatch.getScore()) {
					rank = previousMatch.rank;
				}
			}
			match.rank = rank;
			rank++;
			previousMatch = match;
		}
	}
	
	/**
	 * finds the number of times a certain amino acid is found for each spectrum 
	 * @param matches
	 */
	public static void assignRepeatedPeptideCount(ArrayList<Match> matches) {
		//first error check
		if (matches.size() == 0) return;
		Match.setSortParameter(Match.SORT_BY_SPECTRUM_ID_THEN_PEPTIDE);
		Collections.sort(matches);
		Match match;
		Match previousMatch = matches.get(0);
		int rankCount = 1;
		for (int i = 1; i < matches.size(); i++) {
			//see if these are matches for a different spectrum
			match = matches.get(i);
			if (match.getSpectrum().getId() != previousMatch.getSpectrum().getId()) {
				for (int j = i - rankCount; j < i; j++) {
					matches.get(j).repeatCount = rankCount;
				}
				rankCount = 1;
			} else {
				if (match.getPeptide().equals(previousMatch.getPeptide())) {
					rankCount++;
				} else {
					for (int j = i - rankCount; j < i; j++) {
						matches.get(j).repeatCount = rankCount;
					}
					rankCount = 1;
				}
			}
			previousMatch = match;
		}
	}
	
	/**
	 * NOTE:  this returns early if not DNA search.
	 * 
	 * It removes hits to peptides which are redundant.  For example, a spectrum
	 * has two matches where the peptide is the exact same: same sequence, same start, same direction.
	 * These redundant matches can occasionally come up due to the way large sequences are digested
	 * @param matches
	 */
	public static void removeDuplicateMatches(ArrayList<Match> matches) {
		//first error check
		if (matches.size() == 0) return;
		
		//second error check
		if (!Properties.isSequenceFileDNA) return;
		
		Match.setSortParameter(Match.SORT_BY_SPECTRUM_ID_THEN_PEPTIDE);
		Collections.sort(matches);
		int numberOfMatches = matches.size();
		Match match;
		Match previousMatch = matches.get(0);
		int spectrumID;
		int previousSpectrumID = previousMatch.getSpectrum().getId();
		for (int i = 1; i < numberOfMatches; i++) {
			match = matches.get(i);
			spectrumID = match.getSpectrum().getId();

			if (match.equals(previousMatch) && match.getPeptide().getStartIndex() == previousMatch.getPeptide().getStartIndex() && spectrumID == previousSpectrumID) {
				matches.remove(i);
				i--;
				numberOfMatches--;
			} else {
				previousMatch = match;
				previousSpectrumID = spectrumID;
			}
		}
	}
	

	public static void assignConfidenceValuesToMatches(ArrayList<Match> matches, ArrayList<Spectrum> spectra) {
		/* tally histograms for all spectra so E values can be calculated */
		for (Spectrum spectrum: spectra) {
			spectrum.getEValueCalculator().calculateHistogramProperties();
		}
		
		/* now calculate e values for all matches */
		for (Match match: matches) {
			
			/* determine and set the E, P values */
			match.calculateEValue();
			match.calculatePValue();
			
			/* this is a sanity check for overly confident e values */
			if (match.getEValue() < match.getIMP()) {
				match.setEValue(Double.MAX_VALUE);
			}
		}
	}
	
	/**
	 * Removes matches with poor E values, rank, etc.
	 * @param newMatches
	 * @param matches
	 */
	public static ArrayList<Match> removePoorMatches(ArrayList<Match> matches) {	

		boolean remove;
		ArrayList<Match> paredDownMatches = new ArrayList<Match>();
		for (Match match: matches) {
			
			/* the remove flag initially set to false */
			remove = false;
			
			/* flagging matches with low rank */
			if (match.rank > Properties.maximumNumberOfMatchesForASpectrum) {
				remove = true;
			}
			
			/* flagging matches with high E values */
			if (match.getEValue() > Properties.maxEValue) {
				remove = true;
			}
			
			/* flagging matches with weirdly low E values */
//			if (match.getIMP() > match.getEValue()) {
//				remove = true;
//			}
			
			/* removing the match if it was flagged */
			if (!remove) {
				paredDownMatches.add(match);
			}
			
		}
		
		return paredDownMatches;
	}
	
	private static ArrayList<Peptide> getPeptidesInRegions(ArrayList<Region> regions, ArrayList<Sequence> sequences, boolean isReverse) {
		/* our return */
		ArrayList<Peptide> out = new ArrayList<Peptide>();
		
		/* get peptides along regions */
		for (Sequence sequence: sequences) {
			Sequence_DNA dnaSequence = (Sequence_DNA) sequence;
			for (Region region: regions) {
				if (dnaSequence.equals(region.getSequence())) {
					//TODO get rid of this 500
					out.addAll(dnaSequence.extractPeptidesFromRegion(region.getStartLocation() - 500, region.getStopLocation(), isReverse));
				}
			}
			sequence.reset();
		}
		Collections.sort(out);
		return out;
	}
	
	protected static void printGreeting() {
		U.p("Welcome to Peppy");
		U.p("Protein identification / proteogenomic software.");
		U.p("Developed 2010 by Brian Risk");
		U.p();
		
		/* print some system statistics */
		U.p("number of processors available: " + Runtime.getRuntime().availableProcessors());
		U.p("max available memory: " + (double) memoryUsage.getMax() / (1024 * 1024 * 1024) + " gigabytes");
		U.p();
		
//		U.p("Peppy Copyright (C) 2011 Brian Risk");
//		U.p("This program comes with ABSOLUTELY NO WARRANTY;");

	}
	
	protected static void printFarewell() {
		U.p("max memory used: " + (double) maxMemoryUsed / (1024 * 1024 * 1024) + " gigabytes");
		U.p("done");
	}


	public static ArrayList<Match> reduceMatchesToOnePerSpectrum(ArrayList<Match> matches) {
		Match.setSortParameter(Match.SORT_BY_SPECTRUM_ID_THEN_SCORE);
		Collections.sort(matches);
		int previousID = -1;
		int id;
		ArrayList<Match> out = new ArrayList<Match>(matches.size());
		for (Match match: matches) {
			id = match.getSpectrum().getId();
			if (id != previousID) {
				previousID = id;
				out.add(match);
			}
		}
		return out;
	}
	
	

}
