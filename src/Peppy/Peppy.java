package Peppy;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;

import Reports.HTMLReporter;
import Reports.TextReporter;

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
	
	/* track how much memory we have used */
	static MemoryUsage memoryUsage;
	static long maxMemoryUsed = 0;
	static String propertyFileString = null;
	
	private static final long expiration = 1328292741160L + (1000 * 60 * 60 * 24 * 30);
	
	
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
	
	/**
	 * There may be a multiplicity of sequence and spectral directories.
	 * This iterates through them in every combination.
	 * @param args
	 */
	public static void runPeppy(String [] args) {
		U.startStopwatch();
		
		try {
			
			/* we shall use this variable to track how many searches we performed.
			 * this will be used in the labeling of our report directory
			 */
			int reportIndex = 0;
			
			/* create new report directory */
			String reportDirString;
			if (propertyFileString != null) {
				reportDirString = U.getFileNameWithoutSuffix(new File(propertyFileString));
			} else {
				reportDirString = Properties.spectraDirectoryOrFile.getName() + "_" + System.currentTimeMillis();
			}
			File mainReportDir = new File(Properties.reportDirectory, reportDirString);
			mainReportDir.mkdirs();
			
			/* this will maintain our list of score cutoffs */
			PrintWriter fdrCutoffs = new PrintWriter(new FileWriter (new File(mainReportDir, "FDR cutoffs.txt")));
			
			/* This is what you need to do next:
			 * create an int field inside of peptide
			 * that can act as a general label.  This will let you label peptides
			 * as coming from reverse database, or which search it came from.  In the case of the latter,
			 * that can then be used to keep track of which peptides came where when we lump them all together for
			 * the blind PTM section.
			 */
		
		
			for (File spectraDirectoryOrFile: Properties.spectraDirectoryOrFileList) {
				Properties.spectraDirectoryOrFile = spectraDirectoryOrFile;
				ArrayList<Spectrum> spectra = Spectrum.loadSpectra();
				int spectraSize = spectra.size();
				U.p("loaded " + spectraSize + " spectra.");
				
				ArrayList<Match> allMatchesForSpectralSet = new ArrayList<Match>();
				
				/* iterate through all of our peptide sources */
				for (int sequenceIndex = 0; sequenceIndex < Properties.sequenceDirectoryOrFileList.size(); sequenceIndex++) {
					
					/* set up our sequence data */
					Properties.sequenceDirectoryOrFile = Properties.sequenceDirectoryOrFileList.get(sequenceIndex);
					Properties.isSequenceFileDNA = Properties.isSequenceFileDNAList.get(sequenceIndex);
					ArrayList<Sequence> sequences = Sequence.loadSequenceFiles(Properties.sequenceDirectoryOrFile);
					
					/* increment our report index */
					reportIndex++;
					
					/* create the directory where we will hold this report */
					String reportDirName = reportIndex + " " + spectraDirectoryOrFile.getName() + " - " + Properties.sequenceDirectoryOrFile;
					File reportDir = new File (mainReportDir, reportDirName);
					reportDir.mkdirs();
					
					/* if we are out of spectra (not likely), get out of this loop */
					if (spectraSize == 0) break;
					
					
					/*create an FDR object.  This will define where our score cutoff for this pass is. */
					FDR fdr = new FDR(spectra, reportDir);
					double minimumScore = fdr.getScoreThreshold(Properties.maximumFDR);
					fdrCutoffs.println(reportDirName);
					fdrCutoffs.println(minimumScore);
					fdrCutoffs.println();
					fdr.saveReport();
					
					/* set our minimum score to keep only the confident matches */
					Properties.minimumScore = minimumScore;
					
					/* if we will not find any matches with confidence, skip this round */
					//NOTE:  this will produce no report.  Look out for this when assembling reports!
					if (Properties.minimumScore < 0) continue;
					
					/* get the matches */
					ArrayList<Match> matches = runSearch(spectra, sequences, reportDir);
					
					/* label the matches according to our sequence */
					for (Match match: matches) {
						match.setTrackingIdentifier(sequenceIndex);
					}
					
					/* add these matches to our large pile of all matches for this spectral set */
					allMatchesForSpectralSet.addAll(matches);
					
					/* group spectra that have been identified */
					Hashtable<Integer, Integer> spectrumIDs = new Hashtable<Integer, Integer>(spectra.size());
					for (Match match: matches) {
						spectrumIDs.put(match.getSpectrum().getId(), match.getSpectrum().getId());
					}
					U.p(spectrumIDs.size() + " spectra identified at this step");
					
					/* don't do this if this is the last iteration */
					if (sequenceIndex == Properties.sequenceDirectoryOrFileList.size() - 1) {
						
						/* remove all spectra that appear in our matches */
						Spectrum spectrum;
						for (int spectrumIndex = 0; spectrumIndex < spectra.size(); spectrumIndex++) {
							spectrum = spectra.get(spectrumIndex);
							if (spectrumIDs.get(spectrum.getId()) != null) {
								spectra.remove(spectrumIndex);
								spectrumIndex--;
							}
						}
						double precentReduction =  (1.0 - ((double)spectra.size() / spectraSize));
						U.p("spectra reduced by " + Properties.nfPercent.format(precentReduction) + "%");
						spectraSize = spectra.size();
					}
				}
				
				/* now that we have collected all of our matches from the first pass,
				 * we can collect all the found peptides and perform varible-mod searches on them
				 */
			}
			
			fdrCutoffs.flush();
			fdrCutoffs.close();
		
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		U.stopStopwatch();
	}
	
	
	/**
	 * A "search" is one set of spectra on one peptide database
	 * @param spectra
	 * @param sequences
	 */
	public static ArrayList<Match> runSearch(ArrayList<Spectrum> spectra, ArrayList<Sequence> sequences, File reportDir) {

		//save our properties
		Properties.generatePropertiesFile(reportDir);

		
		/* set up where we will hold all of the matches for our spectra */
		ArrayList<MatchesSpectrum> spectraMatches = new ArrayList<MatchesSpectrum>(spectra.size());
		for (Spectrum spectrum: spectra) {
			spectraMatches.add(new MatchesSpectrum(spectrum));
		}
		
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
			matches = getMatchesWithPeptides(peptides, spectraMatches);
			
		} else {	
			if (Properties.useSequenceRegion) {
				U.p("digesting part of sequence");
				ArrayList<Sequence> oneSequenceList = new ArrayList<Sequence>();
				oneSequenceList.add(sequences.get(0));
				sequences = oneSequenceList;
			}
			
			/* nothing weird.  Just do a normal search */
			matches = getMatches(sequences, spectraMatches);

		}
		
	
//		/* regions (only work with DNA) */
//		if (Properties.isSequenceFileDNA) {	
//			
//			/* rerun the regions analysis */
//			Regions regions = new Regions(matches, sequences);
//			
//			/* creating regions report */
//			regions.createReport(reportDir);
//			
//			/* clear memory */
//			regions.clearRegions();
//		}
		
		U.p("creating text reports");
		TextReporter textReport = new TextReporter(matches, spectra, sequences, reportDir);
		textReport.generateFullReport();
		
		if (Properties.createHTMLReport) {
			U.p("creating HTML reports");
			HTMLReporter report = new HTMLReporter(matches, spectra, sequences, reportDir);
			report.generateFullReport();
		}	
		
		
			
		return matches;
	
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
				propertyFileString = jobFiles.get(i).getName();
				init(args);
				init(jobFiles.get(i).getAbsolutePath());
				try {
					runPeppy(null);
				}
				catch (Exception e) {
					U.p("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
					e.printStackTrace();
					U.p("\r\r");
					U.p("An error occurred. Continuing with the next job...\r\r");
				}
			}
		}
	}


	public static void init(String propertiesFile) {
		/* check to see if this version is too old */
//		if (System.currentTimeMillis() > expiration) {
//			U.p("Welcome to Peppy");
//			U.p("This is an out of date version.  Please access PeppyResearch.com for the current version.");
//			System.exit(0);
//		}
		
		/* this allows us to do our graphics */
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
	
	
	public static ArrayList<Match> getReverseMatches(ArrayList<Sequence> sequences, ArrayList<MatchesSpectrum> spectraMatches) {
		
		/* performing our normal search specifying "true" for isReverse */
		ArrayList<Match> matches = getMatches(sequences, spectraMatches, true);
		
		/* performing the multipass search */
		if (Properties.multipass && Properties.isSequenceFileDNA) {
			matches = multipass(matches, sequences, spectraMatches);
		}
		
		/* label the reverse matches as being from a reverse database */
		for (Match match: matches) {
			match.getPeptide().setDecoy(true);
		}
		
		return matches;
	}
	
	
	/**
	 * Does a normal search, then does multipass
	 * @param sequences
	 * @param spectra
	 * @return
	 */
	public static ArrayList<Match> getMatches(ArrayList<Sequence> sequences, ArrayList<MatchesSpectrum> spectraMatches) {
			
			/* performing our normal search */
			ArrayList<Match> matches = getMatches(sequences, spectraMatches, false);
			
			/* performing the multipass search */
			if (Properties.multipass && Properties.isSequenceFileDNA) {
				matches = multipass(matches, sequences, spectraMatches);
			}
			return matches;
	}

	/**
	 * This is the heart of Peppy where the grand symphony takes place.
	 * 
	 * @param sequences our list of sequences where we will be getting our peptides
	 * @param spectra our list of spectra
	 * @param isReverse if we are doing a normal, forwards search or if this is a null, reverse search
	 * @return
	 */
	private static ArrayList<Match> getMatches(
			ArrayList<Sequence> sequences, 
			ArrayList<MatchesSpectrum> spectraMatches, 
			boolean isReverse) {
			
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
				ArrayList<Peptide> peptides = new ArrayList<Peptide>(Properties.desiredPeptideDatabaseSize);
				peptides.ensureCapacity(Properties.desiredPeptideDatabaseSize);
				
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
				
				/* find the matches for this chunk of peptides */
				(new ScoringThreadServer(peptides, spectraMatches)).findMatches();
				
				/* free up memory */
				peptides.clear();
				System.gc();
					
				/* break if we have covered our last sequence or if we are only using a region */
				if (sequenceIndex == sequences.size() || Properties.useSequenceRegion) {
					break;
				}

			}
			

		return getMatchesFromSpectraMatches(spectraMatches);
	}
	
	
	
	/**
	 * Gets matches where a list of peptides is already derived
	 * @param peptides
	 * @param spectra
	 * @param sequence_DNA
	 * @return
	 */
	public static ArrayList<Match> getMatchesWithPeptides(
			ArrayList<Peptide> peptides, 
			ArrayList<MatchesSpectrum> spectraMatches) {
		
		//This is where the bulk of the processing in long jobs takes
		(new ScoringThreadServer(peptides, spectraMatches)).findMatches();
		ArrayList<Match> matches = getMatchesFromSpectraMatches(spectraMatches);
		
		return matches;
	}
	
	
	public static ArrayList<Match> getMatchesFromSpectraMatches(ArrayList<MatchesSpectrum> spectraMatches) {
		
		/* calculate size of combined matches */
		int size = 0;
		for (MatchesSpectrum spectrumMatches: spectraMatches) {
			size += spectrumMatches.getMatches().size();
		}
		
		/* combine matches into result and return */
		ArrayList<Match> out = new ArrayList<Match>(size);
		for (MatchesSpectrum spectrumMatches: spectraMatches) {
			out.addAll(spectrumMatches.getMatches());
		}
		return out;
	}
	
	
	
	public static ArrayList<Match> reduceMatchesToOnePerSpectrum(ArrayList<Match> matches) {
		Hashtable<String, Match> oneMatchPerSpectrum = new Hashtable<String, Match>(matches.size());
		for (Match match: matches) {
			Match existing = oneMatchPerSpectrum.get(match.getSpectrum().getMD5());
			if (existing == null) {
				oneMatchPerSpectrum.put(match.getSpectrum().getMD5(), match);
			} else {
				if (match.getScore() > existing.getScore()) {
					oneMatchPerSpectrum.put(match.getSpectrum().getMD5(), match);
				}
			}
		}
		ArrayList<Match> out = new ArrayList<Match>(oneMatchPerSpectrum.size());
		Enumeration<Match> e = oneMatchPerSpectrum.elements();
		while (e.hasMoreElements()) out.add(e.nextElement());
		return out;
	}
	


	private static ArrayList<Match> multipass(ArrayList<Match> matches, ArrayList<Sequence> sequences, ArrayList<MatchesSpectrum> spectraMatches) {
		/* regions */
		Regions regions = new Regions(matches, sequences);
		U.p("we found this many regions: " + regions.getRegions().size());
		
		if (regions.getRegions().size() > 0) {
		
			U.p("performing localized PTM search");
			Properties.matchConstructor = new MatchConstructor("Peppy.Match_IMP_VariMod");
			Properties.searchModifications = true;
			
			/* first we need to level the playing field and score the original matches as modification matches */
			//TODO: this can be parallelized
			ArrayList<Match> rescoredMatches = new ArrayList<Match>(matches.size());
			for (Match match: matches) {
				rescoredMatches.add(Properties.matchConstructor.createMatch(match.getSpectrumMatches(), match.getPeptide()));
			}
			matches = rescoredMatches;
			
			/* generating peptides in the regions */
			ArrayList<Peptide> peptidesInRegions = getPeptidesInRegions(regions.getRegions(), sequences, false);
			U.p("there are this many peptides in the regions: " + peptidesInRegions.size());
			
			/* getting matches to the peptides in the regions */
			ArrayList<Match> modMatches = getMatchesWithPeptides(peptidesInRegions, spectraMatches);
			
			/* make room for the modified matches and add them to our main set of matches */
			int matchesNotFoundInUnmodifiedSearchCount = 0;
			for (Match match: modMatches) {
				if (match.isFromModificationSearches()) matchesNotFoundInUnmodifiedSearchCount++;
			}
			matches.ensureCapacity(matches.size() + matchesNotFoundInUnmodifiedSearchCount);
			for (Match match: modMatches) {
				if (match.isFromModificationSearches()) matches.add(match);
			}
			
			/* return things to the way they were */
			Properties.matchConstructor = new MatchConstructor("Peppy.Match_IMP");
			Properties.searchModifications = false;
		}
		
		return matches;
	}


	//TODO this could probably be made much, much faster
	private static ArrayList<Peptide> getPeptidesInRegions(ArrayList<Region> regions, ArrayList<Sequence> sequences, boolean isReverse) {
		/* our return */
		ArrayList<Peptide> out = new ArrayList<Peptide>();
		
		/* determine how many regions we will use to generate our new set of peptides */
		int numberOfRegionsToKeep = Properties.numberOfRegionsToKeep;
		if (numberOfRegionsToKeep > regions.size()) numberOfRegionsToKeep = regions.size();
		
		/* get peptides along regions */
		Region region;
		for (Sequence sequence: sequences) {
			Sequence_DNA dnaSequence = (Sequence_DNA) sequence;
			for (int regionIndex = 0; regionIndex < numberOfRegionsToKeep; regionIndex++) {
				region = regions.get(regionIndex);
				if (dnaSequence.equals(region.getSequence())) {
					
					//TODO get rid of this 500
					out.addAll(dnaSequence.extractPeptidesFromRegion(region.getStartLocation() - 200, region.getStopLocation() + 200, isReverse));
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
		U.p("done");
	}
	
	

}
