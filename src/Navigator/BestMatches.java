package Navigator;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;

import Database.Column;
import Database.Table;
import Peppy.Definitions;
import Peppy.Match_Blank;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.Spectrum;
import Peppy.SpectrumLoader;
import Peppy.U;
import Reports.HTMLPageSpectrum;
import Reports.MatchSVG;
import Reports.UCSC;

/**
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class BestMatches {
	
	/* the sample from which all of the results have been derived */
	private String sampleName;
	
	/* where we keep the best results */
	Hashtable<String, Match> bestMatches = new Hashtable<String, Match>();
	
	/* the results, reduced to the best example for each peptide */
	Hashtable<String, Match> bestPeptides = new Hashtable<String, Match>();
	
	/* every match for each identified peptide sequence */
	Hashtable<String, ArrayList<Match>> allPeptides = new Hashtable<String, ArrayList<Match>>();
	
	/* match types */
	ArrayList<ResultsCategory> resultsCategories = new ArrayList<ResultsCategory>();
	
	/* where the peppy results are stored */
	File resultsFolder;
	
	public static void main(String args[]) {
//		washu();
//		washuChr8();
//		pandey();
//		mayo();
//		ucla();
//		yale();
//		washUPaperOne();
//		washURegionsOfInterest();
//		gm12878();
//		compRefUNC();
//		pccCompRefA();
//		seleno();
//		removedNonsense();
//		carthene();
//		nonsenseSearch();
//		compRefUNCRegionsOfInterest();
//		washUPaperOneRegionAnalysis();
//		yaleEnzymeless();
//		report("/Users/risk2/Documents/workspace/JavaGFS/reports/wei-yang/");
//		moonShoot();
//		multiSetReport("/Users/risk2/Documents/workspace/JavaGFS/reports/lyris-20131201/");
		multiSetReport("/Users/risk2/Documents/workspace/JavaGFS/reports");
//		multiSetReport("/Users/risk2/Documents/workspace/JavaGFS/reports/");
//		multiSetReport("/Users/risk2/Documents/workspace/JavaGFS/reports-complete/MoonShot first pass");
//		createRegionsFile("alex-kentsis", "2 HG19");
//		multiSetReport("/Users/risk2/Documents/workspace/JavaGFS/reports-complete/colorectal - first");
//		multiSetReport("/Users/risk2/Documents/workspace/JavaGFS/reports/scott-nicholson");
//		createRegionsFile("/Users/risk2/PeppyData/CPTAC/reports/", "HG19");
		U.p("done");
	}
	
	public static void multiSetReport(String reportDir) {
		
		ArrayList<File> reportFolders = new ArrayList<File>();
		
		File reports = new File(reportDir);
		File [] allReportFolders = reports.listFiles();
		for (File reportFolder: allReportFolders) {
			if (reportFolder.isDirectory() && !reportFolder.isHidden()) {
				reportFolders.add(reportFolder);
			}
		}
		
		
		/* a list of our BestMatches */
		ArrayList<BestMatches> bestMatches = new ArrayList<BestMatches>();
		
		/* a list of folders to ignore from our results */
		ArrayList<String> direcotryTitlesToIgnore = new ArrayList<String>();
		direcotryTitlesToIgnore.add("varimod");
		direcotryTitlesToIgnore.add("UniProt");
		//direcotryTitlesToIgnore.add("pseudogenes");
		direcotryTitlesToIgnore.add("mouse");
		direcotryTitlesToIgnore.add("UniProt");
		direcotryTitlesToIgnore.add("xeno");
		direcotryTitlesToIgnore.add("germline");
		direcotryTitlesToIgnore.add("gencode");
		direcotryTitlesToIgnore.add("p5");
		direcotryTitlesToIgnore.add("p6");
		
		for (File folder: reportFolders) {
//			BestMatches matches = new BestMatches(folder, ResultsCategory.DNA, direcotryTitlesToIgnore);
			BestMatches matches = new BestMatches(folder, ResultsCategory.PROTEIN, direcotryTitlesToIgnore);
			bestMatches.add(matches);
		}
		
		createUnifiedSamplesReport(bestMatches, "peptideSequence");
	}
	
	
	
	public static void report(String fileName) {
		ArrayList<File> reportFolders = new ArrayList<File>();
		reportFolders.add(new File(fileName));
		
		/* a list of our BestMatches */
		ArrayList<BestMatches> bestMatchesArray = new ArrayList<BestMatches>();
		
		/* a list of folders to ignore from our results */
		ArrayList<String> direcotryTitlesToIgnore = new ArrayList<String>();
		direcotryTitlesToIgnore.add("varimod");
		direcotryTitlesToIgnore.add("UniProt");
		
		for (File folder: reportFolders) {
			BestMatches matches = new BestMatches(folder, -1, direcotryTitlesToIgnore);
			bestMatchesArray.add(matches);
		}
		
		createUnifiedSamplesReport(bestMatchesArray, "peptideSequence");
		U.p("done");
	}
	

	public BestMatches(String sampleName) {
		this.sampleName = sampleName;
	}


	public BestMatches(File resultsFolder, int resultsTypeToAccept, ArrayList<String> direcotryTitlesToIgnore) {
		this( resultsFolder,  resultsTypeToAccept, direcotryTitlesToIgnore, resultsFolder.getName());
		
	}
	 
	public BestMatches(File resultsFolder, int resultsTypeToAccept, ArrayList<String> direcotryTitlesToIgnore, String sampleName) {
		this.sampleName = sampleName;
		this.resultsFolder = resultsFolder;
		
		/* list all the directories in this folder */
		File [] reportFolders = resultsFolder.listFiles();
		U.p(resultsFolder.getName());
		for (File reportFolder: reportFolders) {
			if (!reportFolder.isDirectory()) continue;
			
			/* get out if the report folder is one we are ignoring */
			if (direcotryTitlesToIgnore != null) {
				boolean ignoreDirectory = false;
				for (String directoryTitleToIgnore: direcotryTitlesToIgnore) {
					if (reportFolder.getName().toLowerCase().indexOf(directoryTitleToIgnore.toLowerCase()) != -1) ignoreDirectory = true;
				}
				if (ignoreDirectory) continue;
			}
			
			/* the file where the text report file should be*/
			File textReportFile = new File (reportFolder, "report.txt");
			
			/* get out if there is no report file */
			if (!textReportFile.exists()) continue;
			
			String folderName = reportFolder.getName();
	
			/* extract the database name from the directory name */
			String [] nameComponents = folderName.split("-");
			
			/* database name is the last index */
			String databaseName = nameComponents[nameComponents.length - 1];
			databaseName = databaseName.trim();
			
			/* find out results type (i.e., if this is GENOME or protein) */
			BufferedReader br;
			int resultsType = -1;
			try {
				br = new BufferedReader(new FileReader(textReportFile));
				/* skip the first line */
				br.readLine();
				
				/* > analysis-type: nucleotide */
				String line = br.readLine();
				
				if (line.endsWith("nucleotide")) {resultsType = ResultsCategory.DNA;}
				if (line.endsWith("protein")) {resultsType = ResultsCategory.PROTEIN;}
				
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
			/* get out if we find this is not the correct results type */
			if (resultsTypeToAccept != -1) {
				if (resultsType != resultsTypeToAccept) continue;
			}
	
			/* add these results to the bestMatches */
			ResultsCategory results = new ResultsCategory(databaseName, resultsType);
			results.addFile(textReportFile);
			addMatchType(results);
			
		}	
		process();
		
	}

	
	public static void moonShoot() {
		
		ArrayList<File> reportFolders = new ArrayList<File>();
		
		File reports = new File("reports");
		File [] allReportFolders = reports.listFiles();
		for (File reportFolder: allReportFolders) {
			if (reportFolder.isDirectory()) {
				String [] nameChunks = reportFolder.getName().split("_");
				try {
					Integer.parseInt(nameChunks[0]);
					reportFolders.add(reportFolder);
					
				} catch(NumberFormatException e) {  
			    }
			}
		}
		
		
		/* a list of our BestMatches */
		ArrayList<BestMatches> bestMatches = new ArrayList<BestMatches>();
		
		/* a list of folders to ignore from our results */
		ArrayList<String> direcotryTitlesToIgnore = new ArrayList<String>();
		direcotryTitlesToIgnore.add("varimod");
		direcotryTitlesToIgnore.add("UniProt");
		
		for (File folder: reportFolders) {
//			BestMatches matches = new BestMatches(folder, ResultsCategory.DNA, direcotryTitlesToIgnore);
			BestMatches matches = new BestMatches(folder, ResultsCategory.PROTEIN, direcotryTitlesToIgnore);
			bestMatches.add(matches);
		}
		
		createUnifiedSamplesReport(bestMatches, "peptideSequence");
	}
	
	/**
	 * Used with createRegionsFile to find all folders that match, say HG19
	 */
	public static void findFolderPatterns(File folder, String filePattern, ArrayList<File> filesWithPattern) {
		File [] files = folder.listFiles();
		for (File file: files) {
			if (file.isDirectory()) {
				if (file.getName().indexOf(filePattern) != -1) {
					File reportFile = new File(file, "report.txt");
					if (reportFile.exists()) filesWithPattern.add(reportFile);
				} else {
					findFolderPatterns(file, filePattern, filesWithPattern);
				}
			}
		}
	}
	
	public static void createRegionsFile(String reportsString, String genomeString) {
		
		File reports = new File(reportsString);
		ArrayList<File> filesWithPattern = new ArrayList<File>();
		findFolderPatterns(reports, genomeString, filesWithPattern);

		ArrayList<Match> allMatches = new ArrayList<Match>();
		for (File reportFile: filesWithPattern) {
			allMatches.addAll( Match.loadMatches(reportFile));
		}
		
		
		/* here we are storing all the potential loci
		 * This will be used to create our targeted regions
		 */
		Hashtable<String, Match> targetedRegions = new Hashtable<String, Match>();
		
		for (Match match: allMatches) {
			double startLocus = match.getInt("start");
			//we are rounding this to a thousand so that it gets the general regions
			// thus if three matches are to the same region, it won't produce three regions
			int roundedLocus = (int) (Math.round(startLocus/1000) * 1000);
			String sequenceName = match.getString("sequenceName");
			String targetregion = sequenceName + roundedLocus;
			targetedRegions.put(targetregion, match);
		}
		
		/* creating targeted regions file */
		File targetedRegionsFile = new File("targeted regions file.txt");
		PrintWriter pw;
		try {
			pw = new PrintWriter(new BufferedWriter(new FileWriter(targetedRegionsFile)));
			for(Match match: targetedRegions.values()) {
				int startLocus = match.getInt("start");
				String sequenceName = match.getString("sequenceName");
				if (sequenceName.startsWith(">")) sequenceName = sequenceName.substring(1);
				String peptideSequence = match.getString("peptideSequence");
				pw.println(peptideSequence + "\t" + sequenceName + "\t" + startLocus);
			}
			pw.flush();
			pw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}


	
	
	public static void washUPaperOne() {
		
		ArrayList<File> reportFolders = new ArrayList<File>();
		
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM2-Ellis033/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM2-Ellis041/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM2-Ellis043/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM2-Ellis043/"));
		
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM16-Ellis033/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM16-Ellis041/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM16-Ellis043/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM16-Ellis043/"));
		
		
		/* a list of our BestMatches */
		ArrayList<BestMatches> bestMatches = new ArrayList<BestMatches>();
		
		/* a list of folders to ignore from our results */
		ArrayList<String> direcotryTitlesToIgnore = new ArrayList<String>();
		direcotryTitlesToIgnore.add("varimod");
		
		for (File folder: reportFolders) {
			BestMatches matches = new BestMatches(folder, ResultsCategory.DNA, direcotryTitlesToIgnore);
			bestMatches.add(matches);
		}
		
		createUnifiedSamplesReport(bestMatches, "peptideSequence");
	}
	
	
	
	
	public static void washURegionsOfInterest() {
		
		
		
		ArrayList<File> subtractReportFolders = new ArrayList<File>();
		
		subtractReportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM2-Ellis033/"));
		subtractReportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM2-Ellis041/"));
		subtractReportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM2-Ellis043/"));
		subtractReportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM2-Ellis043/"));
		
		subtractReportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM16-Ellis033/"));
		subtractReportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM16-Ellis041/"));
		subtractReportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM16-Ellis043/"));
		subtractReportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM16-Ellis043/"));
		
//		/*subtract the previous, buggy results */
//		subtractReportFolders.add(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/regions attempt 3/Region-WHIM2-Ellis033/"));
//		subtractReportFolders.add(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/regions attempt 3/Region-WHIM2-Ellis041/"));
//		subtractReportFolders.add(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/regions attempt 3/Region-WHIM2-Ellis043/"));
//		subtractReportFolders.add(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/regions attempt 3/Region-UNC-WHIM2-Ellis043/"));
//		
//		subtractReportFolders.add(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/regions attempt 3/Region-WHIM16-Ellis033/"));
//		subtractReportFolders.add(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/regions attempt 3/Region-WHIM16-Ellis041/"));
//		subtractReportFolders.add(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/regions attempt 3/Region-WHIM16-Ellis043/"));
//		subtractReportFolders.add(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/regions attempt 3/Region-UNC-WHIM16-Ellis043/"));
		
		ArrayList<String> subtractDirecotryTitlesToIgnore = new ArrayList<String>();
		subtractDirecotryTitlesToIgnore.add("varimod");
//		subtractDirecotryTitlesToIgnore.add("mouse");
		
		ArrayList<BestMatches> subtractBestMatchesArray = new ArrayList<BestMatches>();
		for (File folder: subtractReportFolders) {
			BestMatches matches = new BestMatches(folder, -1, subtractDirecotryTitlesToIgnore);
			subtractBestMatchesArray.add(matches);
		}
		
		
		
		
		ArrayList<File> reportFolders = new ArrayList<File>();

		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/RegionsOfInterest/WHIM2-Ellis033/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/RegionsOfInterest/WHIM2-Ellis041/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/RegionsOfInterest/WHIM2-Ellis043/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/RegionsOfInterest/UNC-WHIM2-Ellis043/"));
		
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/RegionsOfInterest/WHIM16-Ellis033/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/RegionsOfInterest/WHIM16-Ellis041/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/RegionsOfInterest/WHIM16-Ellis043/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/RegionsOfInterest/UNC-WHIM16-Ellis043/"));
		
		
		/* a list of our BestMatches */
		ArrayList<BestMatches> bestMatchesArray = new ArrayList<BestMatches>();
		
		/* a list of folders to ignore from our results */
		ArrayList<String> direcotryTitlesToIgnore = new ArrayList<String>();
		direcotryTitlesToIgnore.add("varimod");
		direcotryTitlesToIgnore.add("gencode");
		direcotryTitlesToIgnore.add("subject");
		direcotryTitlesToIgnore.add("xeno");
		direcotryTitlesToIgnore.add("mouse");
		
		for (File folder: reportFolders) {
			BestMatches matches = new BestMatches(folder, -1, direcotryTitlesToIgnore);
			bestMatchesArray.add(matches);
		}
		
		for (BestMatches bmKeep: bestMatchesArray) {
			for (BestMatches bmSubtract: subtractBestMatchesArray) {
				bmKeep.subtractBestMatchesPeptide(bmSubtract);
			}
			U.p(bmKeep.getSampleName() + "\t number of spectra: " + bmKeep.getBestMatches().size());
		}
		
		createUnifiedSamplesReport(bestMatchesArray, "peptideSequence");
	}
	

	
	

	

	
	public static void washu() {
		/* WHIM 2 */
		BestMatches whim2 = new BestMatches("WHIM2");
		U.p("\rloading WHIM2");
		
		/* reference protein */
		/* subject protein */
		/* xeno protein */
		/* contaminant protein */
		/* reference genome */
		/* subject genome */
		/* tumor genome */
		// /Users/risk2/Documents/workspace/JavaGFS/reports/WHIM16-Ellis033/
		
		ResultsCategory whim2ContaminantProtein = new ResultsCategory("ContaminantProtein", ResultsCategory.PROTEIN);
		whim2ContaminantProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/WHIM2/WashU WHIM2 mouse/report.txt"));
		whim2.addMatchType(whim2ContaminantProtein);
		
		/* target protein */
		ResultsCategory whim2ReferenceProtein = new ResultsCategory("ReferenceProtein", ResultsCategory.PROTEIN);
		whim2ReferenceProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/WHIM2/WashU WHIM2 human/report.txt"));
		whim2.addMatchType(whim2ReferenceProtein);
		
		/* reference genome */
		ResultsCategory whim2ReferenceGenome = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
		whim2ReferenceGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/WHIM2/WashU WHIM2 hg19/report.txt"));
		whim2.addMatchType(whim2ReferenceGenome);
		
		/* subject genome */
		ResultsCategory whim2SubjectGenome = new ResultsCategory("SubjectGenome", ResultsCategory.DNA);
		whim2SubjectGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/WHIM2/WashU WHIM2 germline/report.txt"));
		whim2.addMatchType(whim2SubjectGenome);
		
		/* disease genome */
		ResultsCategory whim2DiseaseGenome = new ResultsCategory("DiseaseGenome", ResultsCategory.DNA);
		whim2DiseaseGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/WHIM2/WashU WHIM2 xeno/report.txt"));
		whim2.addMatchType(whim2DiseaseGenome);
		
		
		/* find the best peptides */
		whim2.process();
		whim2.saveReports();
		
		
		/* this list should be ordered by hierarchy of importance; least important first */
		
		/* WHIM 16 */
		BestMatches whim16 = new BestMatches("WHIM16");
		U.p("\rloading WHIM16");
		
		/* contaminant protein */
		ResultsCategory whim16ContaminantProtein = new ResultsCategory("ContaminantProtein", ResultsCategory.PROTEIN);
		whim16ContaminantProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/WHIM16/WashU WHIM16 mouse/report.txt"));
		whim16.addMatchType(whim16ContaminantProtein);
		
		/* target protein */
		ResultsCategory whim16ReferenceProtein = new ResultsCategory("ReferenceProtein", ResultsCategory.PROTEIN);
		whim16ReferenceProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/WHIM16/WashU WHIM16 human/report.txt"));
		whim16.addMatchType(whim16ReferenceProtein);
		
		/* reference genome */
		ResultsCategory whim16ReferenceGenome = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
		whim16ReferenceGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/WHIM16/WashU WHIM16 hg19/report.txt"));
		whim16.addMatchType(whim16ReferenceGenome);
		
		/* subject genome */
		ResultsCategory whim16SubjectGenome = new ResultsCategory("SubjectGenome", ResultsCategory.DNA);
		whim16SubjectGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/WHIM16/WashU WHIM16 germline/report.txt"));
		whim16.addMatchType(whim16SubjectGenome);
		
		/* disease genome */
		ResultsCategory whim16DiseaseGenome = new ResultsCategory("DiseaseGenome", ResultsCategory.DNA);
		whim16DiseaseGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/WHIM16/WashU WHIM16 xeno/report.txt"));
		whim16.addMatchType(whim16DiseaseGenome);
		
		/* find the best peptides */
		whim16.process();
		whim16.saveReports();
		
		
		/* a list of our BestMatches */
		ArrayList<BestMatches> bestMatches = new ArrayList<BestMatches>();
		bestMatches.add(whim2);
		bestMatches.add(whim16);
		
		createUnifiedSamplesReport(bestMatches, "spectrumMD5");
		
		
	}
	
	
	public static void mayo() {
		BestMatches dogan = new BestMatches("Dogan First");
		
		/* Reference protein */
		ResultsCategory doganProtein = new ResultsCategory("ReferenceProtein", ResultsCategory.PROTEIN);
		doganProtein.addFile(new File("/Users/risk2/PeppyData/Mayo/reports/human-mods/report.txt"));
		dogan.addMatchType(doganProtein);
		
		
		/* Reference genome */
		ResultsCategory doganHG19 = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
		doganHG19.addFile(new File("/Users/risk2/PeppyData/Mayo/reports/Mayo dogan hg19/report.txt"));
		dogan.addMatchType(doganHG19);
		

		/* find the best peptides */
		dogan.process();
		dogan.saveReports();

		
		/* a list of our BestMatches */
		ArrayList<BestMatches> bestMatches = new ArrayList<BestMatches>();
		bestMatches.add(dogan);

		createUnifiedSamplesReport(bestMatches, "peptideSequence");
	}
	
	
	
	public static void ucla() {
		/* Creinhardtii_169 */
		BestMatches controlSet = new BestMatches("ucla loo control");
		BestMatches deprivedSet = new BestMatches("ucla loo deprived");
//		BestMatches syntrofSet = new BestMatches("ucla loo syntrof");
		
		/* control */
		
		/* contaminant protein */
		ResultsCategory controlContaminants = new ResultsCategory("ContaminantProtein", ResultsCategory.PROTEIN);
		controlContaminants.addFile(new File("/Users/risk2/PeppyData/ucla/reports/ucla-Chlamy control-contaminants/report.txt"));
		controlSet.addMatchType(controlContaminants);
		
//		/* reference protein */
//		ResultsCategory controlReferenceProtein = new ResultsCategory("ReferenceProtein", ResultsCategory.PROTEIN);
//		controlReferenceProtein.addFile(new File("/Users/risk2/PeppyData/ucla/reports/ucla-Chlamy control-protein/report.txt"));
//		controlSet.addMatchType(controlReferenceProtein);
		
		/* reference genome */
		ResultsCategory controlReferenceGenome = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
		controlReferenceGenome.addFile(new File("/Users/risk2/PeppyData/ucla/reports/ucla-Chlamy control-dna/report.txt"));
		controlSet.addMatchType(controlReferenceGenome);
		
		
		
		/* nitrogen deprived */
		
		/* contaminant protein */
		ResultsCategory deprivedContaminants = new ResultsCategory("ContaminantProtein", ResultsCategory.PROTEIN);
		deprivedContaminants.addFile(new File("/Users/risk2/PeppyData/ucla/reports/ucla-Chlamy nitrogen depletion-contaminants/report.txt"));
		deprivedSet.addMatchType(deprivedContaminants);
		
//		/* reference protein */
//		ResultsCategory deprivedReferenceProtein = new ResultsCategory("ReferenceProtein", ResultsCategory.PROTEIN);
//		deprivedReferenceProtein.addFile(new File("/Users/risk2/PeppyData/ucla/reports/ucla-Chlamy nitrogen depletion-protein/report.txt"));
//		deprivedSet.addMatchType(deprivedReferenceProtein);
		
		/* reference genome */
		ResultsCategory deprivedReferenceGenome = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
		deprivedReferenceGenome.addFile(new File("/Users/risk2/PeppyData/ucla/reports/ucla-Chlamy nitrogen depletion-dna/report.txt"));
		deprivedSet.addMatchType(deprivedReferenceGenome);
		
		
		
		
//		/* syntrof */
//		
//		/* contaminant protein */
//		ResultsCategory syntrofContaminants = new ResultsCategory("ContaminantProtein", ResultsCategory.PROTEIN);
//		syntrofContaminants.addFile(new File("/Users/risk2/PeppyData/ucla/reports/ucla-syntrof-contaminants/report.txt"));
//		syntrofSet.addMatchType(syntrofContaminants);
//		
//		/* reference protein */
//		ResultsCategory syntrofReferenceProtein = new ResultsCategory("ReferenceProtein", ResultsCategory.PROTEIN);
//		syntrofReferenceProtein.addFile(new File("/Users/risk2/PeppyData/ucla/reports/ucla-syntrof-protein/report.txt"));
//		syntrofSet.addMatchType(syntrofReferenceProtein);
//		
//		/* reference genome */
//		ResultsCategory syntrofReferenceGenome = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
//		syntrofReferenceGenome.addFile(new File("/Users/risk2/PeppyData/ucla/reports/ucla-syntrof-dna/report.txt"));
//		syntrofSet.addMatchType(syntrofReferenceGenome);
		
		
		

		/* find the best peptides */
		controlSet.process();
		controlSet.saveReports();
		deprivedSet.process();
		deprivedSet.saveReports();
//		syntrofSet.process();
//		syntrofSet.saveReports();

		
		/* a list of our BestMatches */
		ArrayList<BestMatches> bestMatches = new ArrayList<BestMatches>();
		bestMatches.add(controlSet);
		bestMatches.add(deprivedSet);
//		bestMatches.add(syntrofSet);
		
		createUnifiedSamplesReport(bestMatches, "peptideSequence");
	}
	

	
	
	
	public static void pccCompRefA() {
		ArrayList<File> reportFolders = new ArrayList<File>();

		
//		reportFolders.add(new File("reports/BI-P5"));
//		reportFolders.add(new File("reports/JHU-P5"));
//		reportFolders.add(new File("reports/VAND-P5"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM2-Ellis033/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/WashU_045_046_P5"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM2-Ellis041/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM2-Ellis043/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM2-Ellis043/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/CompRef_Proteome_VU_B1_P5"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/CompRef_Proteome_VU_B2_P5"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/CompRef_Proteome_VU_B3_P5"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/CompRef_Proteome_JHUC_P5AB"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM2-CompRef/"));
		
//		reportFolders.add(new File("reports/JHU-P6"));
//		reportFolders.add(new File("reports/VAND-P6"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM16-Ellis033/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/WashU_045_046_P6"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM16-Ellis041/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM16-Ellis043/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM16-Ellis043/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM16-CompRef/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/CompRef_Proteome_VU_A2_P6"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/CompRef_Proteome_VU_A3_P6"));		
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/CompRef_Proteome_JHUC_P6ST"));
		
		
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/CompRef_Proteome_PNNL"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/CompRef_Proteome_JHUC_iTRAQ"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/CompRef_Proteome_BI"));
//		reportFolders.add(new File("reports/PNNL"));
//		reportFolders.add(new File("reports/broad"));

		
		
		
		/* a list of our BestMatches */
		ArrayList<BestMatches> bestMatchesArray = new ArrayList<BestMatches>();
		
		/* a list of folders to ignore from our results */
		ArrayList<String> direcotryTitlesToIgnore = new ArrayList<String>();
		direcotryTitlesToIgnore.add("varimod");
//		direcotryTitlesToIgnore.add("HG19");
		direcotryTitlesToIgnore.add("mouse");
		direcotryTitlesToIgnore.add("xeno");
		direcotryTitlesToIgnore.add("gencode");
		direcotryTitlesToIgnore.add("personal");
		direcotryTitlesToIgnore.add("germline");
		direcotryTitlesToIgnore.add("subject");
		
		for (File folder: reportFolders) {
			BestMatches matches = new BestMatches(folder, -1, direcotryTitlesToIgnore);
			bestMatchesArray.add(matches);
		}
		
		createUnifiedSamplesReport(bestMatchesArray, "peptideSequence");
	}
	

	
	
	

	
	
	
	
	

	

	
	
	public static void createUnifiedSamplesReport(ArrayList<BestMatches> bestMatches, String keyName) {
		File unifiedReportsDir = new File("unifiedReports");
		unifiedReportsDir.mkdir();
		createUnifiedSamplesReport(bestMatches, keyName, unifiedReportsDir);
	}
	
	public static void createUnifiedSamplesReport(ArrayList<BestMatches> bestMatches, String keyName, File parentFolder) {
		boolean makeSpectrumPage = true;
		boolean makeGeneReport = true;
		boolean fullColumns = false;
		boolean showUCSC = true;
		
		/* label the matches for their respective samples */
		Enumeration<Match> elements;
		for (BestMatches bm: bestMatches) {
			elements = bm.getBestMatches().elements();
			while (elements.hasMoreElements()) elements.nextElement().set("sampleName", bm.getSampleName());
		}

		
		/* create tables of results for merging */
		Column key = new Column(keyName, String.class);
		ArrayList<Table> tables = new ArrayList<Table>(bestMatches.size());
		for (BestMatches bm: bestMatches) {
			Table table = new Table(bm.getSampleName(), key);
			elements = bm.getBestMatches().elements();
			while (elements.hasMoreElements()) table.add(elements.nextElement());
			tables.add(table);
		}
		
		/* merge the results into one table */
		Table merged = tables.get(0);
		for (int i = 1; i < tables.size(); i++) {
			Table subMerge = merged.combine(tables.get(i));
			merged = subMerge;
		}
		
		ArrayList<String> allKeys = merged.getKeyValues();

		/* get list of intersecting peptides */
		ArrayList<MergeMatchHolder> pmmhs = new ArrayList<MergeMatchHolder>(); 
		for (String keyValue: allKeys) {
			MergeMatchHolder holder = new MergeMatchHolder(keyName);
			for (BestMatches bm: bestMatches) {
				
				/* change this to getBestMatches when doing spectrumMD5 */
				Match match = bm.getBestPeptideMatches().get(keyValue);
				if (match != null) {
					holder.put(bm.getSampleName(), match);
				}
			}
			pmmhs.add(holder);
		}
		Collections.sort(pmmhs);
		
		for (Table table: tables) {
			U.p(table.getName() + ": " + table.size());
		}
		U.p("merged: " + merged.size());
		
		ArrayList<Region> geneRegions = new ArrayList<Region>();
		if (makeGeneReport) {
			U.p("loading gencode transcripts");
			try {
				File gencodeFile = new File("resources/gencode/gencodeReduced.gtf");
				BufferedReader gencodeReader = new BufferedReader(new FileReader(gencodeFile));
				
				
				String line = gencodeReader.readLine();
				while (line != null) {
					String [] chunks1 = line.split("\t");
					String [] chunks2 = chunks1[8].split(";");
					Region region = new Region();
					region.setSequence(chunks1[0]);
					region.setStart(Integer.parseInt(chunks1[3]));
					region.setStop(Integer.parseInt(chunks1[4]));
					if (chunks1[6].equals("-")) region.setForwards(false);
					String geneName = chunks2[4].substring(12, chunks2[4].length() - 1);
					region.setName(geneName);

					String transcriptType = chunks2[5].substring(18, chunks2[5].length() - 1);
					region.setDescription(transcriptType);

					
					geneRegions.add(region);
					line = gencodeReader.readLine();
				}
				
				gencodeReader.close();

				
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		
		/* write merged table */
		try {
			
			/* create a BED file as well */
			PrintWriter bedWriter = new PrintWriter(new BufferedWriter(new FileWriter(new File(parentFolder, "unified_bed.txt"))));
			bedWriter.println("track name=\"Peppy results\" description=\"Peppy results\" visibility=2 itemRgb=\"On\"");
			
			/* header */
			PrintWriter matchWriter = new PrintWriter(new BufferedWriter(new FileWriter(new File(parentFolder, "unified_report.html"))));
			matchWriter.println("<html>");
			matchWriter.println("<head>");
			matchWriter.println("" +
					"<style type=\"text/css\">" +
					"<!--" +
					"@import url(\"http://proteomics.me/resources/reports/style.css\");" +
					"@import url(\"http://proteomics.me/resources/reports/sortable.css\");" +
					"-->" +
					"</style>" +
					"<script type=\"text/javascript\" src=\"http://proteomics.me/resources/reports/sortable.js\"></script>"
					);
			matchWriter.println("</head>");
			matchWriter.println("<body>");
			
			
			
			/* the table */
			matchWriter.println("<table class=\"sortable\" id=\"box-table-a\" >");
			matchWriter.println("<tr>");
			matchWriter.println("<th>#</th>");
			matchWriter.println("<th>peptide</th>");
//			matchWriter.println("<th>*</th>");
			if (makeGeneReport) {
				matchWriter.println("<th>Gene</th>");
//				matchWriter.println("<th>GeneType</th>");
				matchWriter.println("<th>onco</th>");
			}
			
			if (fullColumns) {
				for (int i = 0; i < bestMatches.size(); i++) {
					matchWriter.println("<th>" + (i + 1) + "</th>");
				}
			}
			matchWriter.println("<th>set count</th>");
			matchWriter.println("<th>spectrum count</th>");
//			matchWriter.println("<th>score total</th>");
//			matchWriter.println("<th>P5</th>");
//			matchWriter.println("<th>P6</th>");
//			matchWriter.println("<th>P5/P6</th>");
//			matchWriter.println("<th>single WHIM</th>");
//			matchWriter.println("<th>NIST</th>");
			if (showUCSC) matchWriter.println("<th>UCSC</th>");
			matchWriter.println("<th>unique location</th>");
			matchWriter.println("<th>sequence</th>");
			matchWriter.println("<th>start</th>");
			matchWriter.println("<th>stop</th>");
			matchWriter.println("<th>strand</th>");
			matchWriter.println("<th>notes</th>");
			
			matchWriter.println("</tr>");
			
			int counter = 0;
			for (MergeMatchHolder holder: pmmhs) {
				counter++;
				
				Match anyMatch = holder.get();
				if (anyMatch == null) continue;
				
				/*
				 * HACK
				 * Skipping OPLAH matches or other known proteins
				 */
				int start = anyMatch.getInt("start") ;
				String peptideSequence = anyMatch.getString("peptideSequence");
				if (anyMatch.getString("sequenceName").equals(">chr8")) {
					if (start > 145106167 && start < 145115584) {
						continue;
					}
				}
				if (start == 77771407) continue;
				if (start == 143767494) continue;
				if (start == 73959542) continue;
				if (start == 8926467) continue;
				if (start == 69349598) continue;
				if (start == 77711364) continue;
				if (start == 118948882) continue;
				if (start == 91033136) continue;
				if (start == 7479940) continue;
				if (start == 104335609) continue;
				if (start == 57108181) continue;
				if (start == 1571780) continue;
				if (start == 1571783) continue;
				if (start == 94800588) continue;
				if (start == 94800941) continue;
				if (start == 73958003) continue;
				if (start == 73958477) continue;
				if (start == 73958561) continue;
				if (start == 73959278) continue;
				if (start == 73959308) continue;
				if (start == 73959470) continue;
				if (start == 74689310) continue;
				if (start == 53321621) continue;
				if (start == 53326161) continue;
				if (start == 53326290) continue;
				if (start == 94294661) continue;
				if (start == 94294859) continue;
				if (start == 64099187) continue;
				if (start == 64099649) continue;
				if (start == 57038572) continue;
				if (start == 57033784) continue;
				if (start == 77771407) continue;
				if (start == 77771407) continue;
				if (start == 77771407) continue;
				if (start == 77771407) continue;
				if (start == 77771407) continue;
				if (start == 77771407) continue;
				if (start == 77771407) continue;
				if (start == 77771407) continue;
				if (start == 77771407) continue;
				if (start == 77771407) continue;
				if (start == 77771407) continue;
				if (peptideSequence.equals("ALAEGPGAEGPR")) continue;
				if (peptideSequence.equals("SMFLQMGTTAGPSGEAGGGAGLAEAPR")) continue;
				if (peptideSequence.equals("ASQNSFRIEYDTFGELK")) continue;
				if (peptideSequence.equals("ALQEASEAYLVGLFEDTNLCAIHAK")) continue;
				if (peptideSequence.equals("AGLGQPLPQEELDAMIR")) continue;
				if (peptideSequence.equals("MCNEFFEGFPDK")) continue;
				if (peptideSequence.equals("RSDTMAGGGGSSDGSGR")) continue;
				if (peptideSequence.equals("VSLGGFEITPPVVLR")) continue;
				if (peptideSequence.equals("WGAATTSPAASDAR")) continue;
				if (peptideSequence.equals("MVTGDHPITAK")) continue;
				if (peptideSequence.equals("MANIQPQMLVAGATSIAR")) continue;
				if (peptideSequence.equals("TALLDAAGVASLLTTAEVVVTEIPK")) continue;
				if (peptideSequence.equals("MMDLQHGSLFLR")) continue;
				if (peptideSequence.equals("MDSGDGVTHTVPIYEGYALPHAILR")) continue;
				if (peptideSequence.equals("SLASGLVPAAPPK")) continue;
				if (peptideSequence.equals("MCTCPSGQIAPSCGSR")) continue;
				if (peptideSequence.equals("AVTFIGNSTAIQELFK")) continue;
				if (peptideSequence.equals("MVEPIQGEAGVVVPDPGYLMGVR")) continue;
				if (peptideSequence.equals("LLTPGACSSEVPSAVPSR")) continue;
				if (peptideSequence.equals("IAQTVTAVAGLTSYPFDTVR")) continue;
				if (peptideSequence.equals("MVTPGHACTQK")) continue;
				if (peptideSequence.equals("ALDIEIATYRK")) continue;
				if (peptideSequence.equals("PVDFTGYWK")) continue;
				if (peptideSequence.equals("TYELLNCDK")) continue;
				if (peptideSequence.equals("SSTDSLPGELR")) continue;
				if (peptideSequence.equals("LPGNISSLNVECR")) continue;
				if (peptideSequence.equals("GETLAETWTR")) continue;
				if (peptideSequence.equals("GPGAPAASSPTQK")) continue;
				if (peptideSequence.equals("HGEQDPSLWEQALSYFAR")) continue;
				if (peptideSequence.equals("GESDDSILR")) continue;
				if (peptideSequence.equals("SGFFLFCSEFRPK")) continue;
				if (peptideSequence.equals("TALLDAAGVASLLTTAEVVVTEIPKEEK")) continue;
				if (peptideSequence.equals("VTDYFVVAGLTDTSTLLDQEINRLDTK")) continue;
				if (peptideSequence.equals("APITDIAIIIK")) continue;
				if (peptideSequence.equals("LIPGCEVILATPYGR")) continue;
				if (peptideSequence.equals("AVTLSALIYCLR")) continue;
				if (peptideSequence.equals("ALEAVATEVNSFLTNGPCPASPLSLEEVAMGFVR")) continue;
				if (peptideSequence.equals("XXXXXX")) continue;
				if (peptideSequence.equals("XXXXXX")) continue;
				if (peptideSequence.equals("XXXXXX")) continue;
				if (peptideSequence.equals("XXXXXX")) continue;
				if (peptideSequence.equals("XXXXXX")) continue;
				if (peptideSequence.equals("XXXXXX")) continue;
				
				
				//MD Anderson false positives
				if (peptideSequence.equals("NDDIPEQDSLGLSNLQK")) continue; //NIST
				if (peptideSequence.equals("GIHETTFNSIMK")) continue; //NIST
				if (peptideSequence.equals("TEQGPTGVTMTSNPITWGQIK")) continue; //NIST
				if (peptideSequence.equals("LETEIEALKEELLFMK")) continue; //NIST
				if (peptideSequence.equals("MTTVHAITATQK")) continue; //NIST
				if (peptideSequence.equals("LVGLFEDTNLCAIHAK")) continue; //NIST
				if (peptideSequence.equals("AGPNTNGSQFFICTAK")) continue;
				if (peptideSequence.equals("EEEEEEEEDQEGDR")) continue;
				if (peptideSequence.equals("PSAAAAAAAAGR")) continue; //100% blast
				if (peptideSequence.equals("SCGIHETTFNSIMK")) continue; //100% blast
				if (peptideSequence.equals("LTADLSAETLQAR")) continue; //NIST
				if (peptideSequence.equals("GAYGGGYGGYDDYGGYNDGYGFGSDR")) continue; //100% blast
				if (peptideSequence.equals("ILELAGNAAR")) continue; //NIST
				if (peptideSequence.equals("PLADAHLTK")) continue; //NIST
				if (peptideSequence.equals("FIAVGYVDDTEFVR")) continue; //NIST
				if (peptideSequence.equals("MLNNPLFAGNPQLQEQMR")) continue; //100% blast
				if (peptideSequence.equals("LENTLQDIIQENFPNLAR")) continue; //NIST
				if (peptideSequence.equals("QLMSYPLSEGQLDQK")) continue; //100% blast
				if (peptideSequence.equals("LCYVALDFEQEMATAASSSSLEK")) continue; //NIST
				if (peptideSequence.equals("EAYLVGLFEDTNLCAIHAK")) continue; //NIST
				if (peptideSequence.equals("LPDGQVITIGNER	")) continue; //NIST
				if (peptideSequence.equals("VYAGAGGSGSR")) continue; //semi-tryp of krt18
				if (peptideSequence.equals("LPDGQVITIGNER")) continue; //NIST semi
				if (peptideSequence.equals("VQVEGYNYTGMGNSTNK")) continue; //UCSC semi tryp
				if (peptideSequence.equals("FCGMFCPVEGSSENK")) continue; //fully tryp
				if (peptideSequence.equals("TIDFDSLSVGR")) continue; // fully tryp
				if (peptideSequence.equals("VDPNSPVLLEDPVLCALAK")) continue; //UCSC semi tryp
				if (peptideSequence.equals("VFEFQLTSEEMK")) continue;  //UCSC semi tryp
				if (peptideSequence.equals("IQAVLLPK")) continue; //blast 100%
				if (peptideSequence.equals("DDCGDGSDEASCPK")) continue; //UCSC fully tryp
				if (peptideSequence.equals("KVEEEDEEEEEEEEEEEEEEDE")) continue; //UCSC fully tryp
				if (peptideSequence.equals("ENEKDLEEAEEYK")) continue; //semi tryptic
				if (peptideSequence.equals("NPVDILTYVAWK")) continue; //100% blast
				if (peptideSequence.equals("YLECISCGSSDMSCER")) continue; //UCSC fully tryptic
				if (peptideSequence.equals("SASSGAEGDVSSEREP")) continue; //UCSC fully tryptic
				if (peptideSequence.equals("AAEEAFVNDIDESSPGTEWER")) continue; //UCSC fully tryp
				if (peptideSequence.equals("EQLNGILLHLESELAQTR")) continue; //blast 100%
				if (peptideSequence.equals("SYSQAGITETEWTSGSSK")) continue;//100% UCSC
				if (peptideSequence.equals("GLMTTVHAITATQK")) continue; //100% blast
				if (peptideSequence.equals("XXXXXX")) continue;
				if (peptideSequence.equals("XXXXXX")) continue;
				if (peptideSequence.equals("XXXXXX")) continue;
				if (peptideSequence.equals("XXXXXX")) continue;
				if (peptideSequence.equals("XXXXXX")) continue;
				if (peptideSequence.equals("XXXXXX")) continue;
				if (peptideSequence.equals("XXXXXX")) continue;
				if (peptideSequence.equals("XXXXXX")) continue;
				if (peptideSequence.equals("XXXXXX")) continue;
				if (peptideSequence.equals("XXXXXX")) continue;
				
				
				//blast distance column
				
				
				
				//MD Anderson notes
				String note = "";
				if (peptideSequence.equals("TEQGPTGVTMTSNPITWGQIK")) note = "very unique sequence";
				if (peptideSequence.equals("NISNASCTTNCIAPLAK")) note = "deamidated, semitryptic IISNAS... from GAPDH";
				if (peptideSequence.equals("VDNDEDEHQLSLR")) note = "deamidation";
				if (peptideSequence.equals("DFIMDNCEELIPEYLNFIR")) note = "canonical is VDFIMDN...fully, tryp.  V+Ox = D";
				if (peptideSequence.equals("MIPILQQIIFR")) note = "blast somewhat similar to Filensin, lense protein";
				if (peptideSequence.equals("LLEDGEDFNLGDALDSSN")) note = "c-term semitryptic, krt18";
				if (peptideSequence.equals("WDYPEGTPNGGSTTLPSAPPPASAGLK")) note = "fairly unique sequene; LPSAPPPA matches atrophin";
				if (peptideSequence.equals("RLLEDGEDFNLGDALDSSN")) note = "100%, but fully non-tryptic, krt18";//keratin 18: http://www.ncbi.nlm.nih.gov/protein/4557888?report=genbank&log$=protalign&blast_rank=2&RID=8FG36HD101R
				if (peptideSequence.equals("EAAAAAAAAHQQALR")) note = "predicted"; //http://www.ncbi.nlm.nih.gov/protein/530404911?report=genbank&log$=protalign&blast_rank=1&RID=8HN9XYK201R
				if (peptideSequence.equals("MVVVWSTGEESLR")) note = "cDNA"; //	http://peptide.nist.gov/browser/peptide_stat.php?description=IT&organism=human&pep_seq=MVVVWSTGEESLR
				if (peptideSequence.equals("IAVPIILR")) note = "BLAST 2";
				if (peptideSequence.equals("NSTLSEPGSGR")) note = "L could be P+Oxy";
				if (peptideSequence.equals("GGQLHAAGGAGR")) note = "frameshift";
				if (peptideSequence.equals("AAAAAAAAAAAAAAAK")) note = "ASAAAAAAAAAAAAAK, -16";
				if (peptideSequence.equals("EAGAGAEAAAGSAR")) note = "BLAST 3";
				if (peptideSequence.equals("WTQLFILER")) note = "BLAST 2";
				if (peptideSequence.equals("IDAIIDLSK")) note = "BLAST 1, but weird difference.  primary sequence, SYSQAGITETEWTSGSSK, is known";
				if (peptideSequence.equals("EANIIFLK")) note = "BLAST 2";
				if (peptideSequence.equals("GDLSIAAIK")) note = "BLAST 2";
				if (peptideSequence.equals("XXXXXX")) note = "xxxxxxxxxx";
				if (peptideSequence.equals("XXXXXX")) note = "xxxxxxxxxx";
				if (peptideSequence.equals("XXXXXX")) note = "xxxxxxxxxx";
				if (peptideSequence.equals("XXXXXX")) note = "xxxxxxxxxx";
				if (peptideSequence.equals("XXXXXX")) note = "xxxxxxxxxx";
				if (peptideSequence.equals("XXXXXX")) note = "xxxxxxxxxx";
				if (peptideSequence.equals("XXXXXX")) note = "xxxxxxxxxx";
				if (peptideSequence.equals("XXXXXX")) note = "xxxxxxxxxx";

				
			
				
				if (peptideSequence.equals("AGGDLPLQPQPGGAAAR")) note = "frame shift";
				if (peptideSequence.equals("SGLAACTAPLPGSR")) note = "reverse frame";
				if (peptideSequence.equals("DAEQEEEVQR")) note = "frame shift";
				if (peptideSequence.equals("ATPGHTGCLSPGCPDQPAR")) note = "frame shift";
				if (peptideSequence.equals("AQIASSGLTVEVDAPK")) note = "mutated stop?";
				if (peptideSequence.equals("XXXXXX")) note = "xxxxxxxxxx";
				if (peptideSequence.equals("XXXXXX")) note = "xxxxxxxxxx";
				if (peptideSequence.equals("XXXXXX")) note = "xxxxxxxxxx";

				
				
				
				
				
				
				
				
				/*
				 * Marking matches being reported in methods paper
				 */
				boolean earlyEllis = false;
				if (peptideSequence.equals("DLADERALVDVIEDK")) earlyEllis = true;
				if (peptideSequence.equals("EQNEAEVASLNR")) earlyEllis = true;
				if (peptideSequence.equals("GEMMDRQHGSLFLR")) earlyEllis = true;
				if (peptideSequence.equals("IGGIASVPVGR")) earlyEllis = true;
				if (peptideSequence.equals("IISNASCTTNCLAPRAK")) earlyEllis = true;
				if (peptideSequence.equals("IMNMIR")) earlyEllis = true;
				if (peptideSequence.equals("VETGFLKPGMVVTVAPVNVTTEVK")) earlyEllis = true;
				if (peptideSequence.equals("ALEMEDSQLCK")) earlyEllis = true;
				if (peptideSequence.equals("LIDEVLK")) earlyEllis = true;
				if (peptideSequence.equals("ARALEMEDSQLCK")) earlyEllis = true;
				if (peptideSequence.equals("ASGVSAAAPGER")) earlyEllis = true;
				if (peptideSequence.equals("DLADELALVDVNEDK")) earlyEllis = true;
				if (peptideSequence.equals("DLADELDLVDVIEDK")) earlyEllis = true;
				if (peptideSequence.equals("EMMDLHHGSLFLR")) earlyEllis = true;
				if (peptideSequence.equals("GEMMDLEHGSLFLR")) earlyEllis = true;
				if (peptideSequence.equals("GLDSLRPLR")) earlyEllis = true;
				if (peptideSequence.equals("GQMMDLQHGSLFLR")) earlyEllis = true;
				if (peptideSequence.equals("GWMMDLQHGSLFLR")) earlyEllis = true;
				if (peptideSequence.equals("HIAEEGDR")) earlyEllis = true;
				if (peptideSequence.equals("ITVVGVGAVGMACANSILMK")) earlyEllis = true;
				if (peptideSequence.equals("LKGEMMDLEHGSLFLR")) earlyEllis = true;
				if (peptideSequence.equals("LLIVSDPVDILTYVAWK")) earlyEllis = true;
				if (peptideSequence.equals("LLIVSNPVDILTYMAWK")) earlyEllis = true;
				if (peptideSequence.equals("LLIVSNPVDILTYVDWK")) earlyEllis = true;
				if (peptideSequence.equals("LLIVSNPVDNLTYVAWK")) earlyEllis = true;
				if (peptideSequence.equals("LLQSIAVK")) earlyEllis = true;
				if (peptideSequence.equals("LPGAAGNR")) earlyEllis = true;
				if (peptideSequence.equals("LQGEMMDLQHGSLFLRTPK")) earlyEllis = true;
				if (peptideSequence.equals("LRIVSNPVDILTYVAWK")) earlyEllis = true;
				if (peptideSequence.equals("NHVIGSVCNLDSAR")) earlyEllis = true;
				if (peptideSequence.equals("QGISLLLR")) earlyEllis = true;
				if (peptideSequence.equals("SMVASGSELGK")) earlyEllis = true;
				if (peptideSequence.equals("SYHEEFNPPQEPMKDDITGEPLIR")) earlyEllis = true;
				if (peptideSequence.equals("THINIIVTGHVDSGK")) earlyEllis = true;
				if (peptideSequence.equals("TQTAAAAAAGGVGGGGGAMGGLASGGDVEPGLPVEVR")) earlyEllis = true;
				if (peptideSequence.equals("VIGSGSNLDSAR")) earlyEllis = true;
				if (peptideSequence.equals("AGPNTNGSQFFICTAK")) earlyEllis = true;
				if (peptideSequence.equals("GIMNSFVNDIFER")) earlyEllis = true;
				if (peptideSequence.equals("MGIMNSFVNDIFER")) earlyEllis = true;
				if (peptideSequence.equals("VPCILGQNGISDLVK")) earlyEllis = true;
				if (peptideSequence.equals("ALEMENSQLCK")) earlyEllis = true;
				if (peptideSequence.equals("APPATNPLATMSVELEEALPVTTAEGMAK")) earlyEllis = true;
				if (peptideSequence.equals("AVAGASAMFAGLQDLGVANGEDLK")) earlyEllis = true;
				if (peptideSequence.equals("DGGLSIPHSTK")) earlyEllis = true;
				if (peptideSequence.equals("DIPHMDIEALK")) earlyEllis = true;
				if (peptideSequence.equals("DRLPLAVVGSDIITEVNGK")) earlyEllis = true;
				if (peptideSequence.equals("GLMTTVHAITATQK")) earlyEllis = true;
				if (peptideSequence.equals("HTGPGILSMANAGLNTNGSQFFICTAK")) earlyEllis = true;
				if (peptideSequence.equals("HTGPGILSMANAGSNINGSQFFICTAK")) earlyEllis = true;
				if (peptideSequence.equals("IEDGNDFGVAIQDK")) earlyEllis = true;
				if (peptideSequence.equals("KIGYNPDTVAFVPISGWNDDNMLEPSANMPWFK")) earlyEllis = true;
				if (peptideSequence.equals("KIGYNPDTVAFVPISGWNGNNMLEPSANMPWFK")) earlyEllis = true;
				if (peptideSequence.equals("KITIANCGQLE")) earlyEllis = true;
				if (peptideSequence.equals("KLGEMWNDTAADDK")) earlyEllis = true;
				if (peptideSequence.equals("KLLMMAGINDCYTSAR")) earlyEllis = true;
				if (peptideSequence.equals("LETEIEALKEELLFMK")) earlyEllis = true;
				if (peptideSequence.equals("LSAIYGGTYMLNKPIEEIIVQDGK")) earlyEllis = true;
				if (peptideSequence.equals("MAGIDDCYTSAR")) earlyEllis = true;
				if (peptideSequence.equals("MELHGEGSSSGK")) earlyEllis = true;
				if (peptideSequence.equals("MGLYHGQVLCK")) earlyEllis = true;
				if (peptideSequence.equals("MTTVHAITATQK")) earlyEllis = true;
				if (peptideSequence.equals("NIELFLNSKIK")) earlyEllis = true;
				if (peptideSequence.equals("NLDVGANIFNGNLDPEIDEK")) earlyEllis = true;
				if (peptideSequence.equals("NPALVELLK")) earlyEllis = true;
				if (peptideSequence.equals("PFIDLNYMVYMFEYDSTHGK")) earlyEllis = true;
				if (peptideSequence.equals("PIEEIIVQDGK")) earlyEllis = true;
				if (peptideSequence.equals("QTTHFVEGGDAGNREDQINR")) earlyEllis = true;
				if (peptideSequence.equals("RSDTMAGGGGSSDGSGR")) earlyEllis = true;
				if (peptideSequence.equals("RSELEAALQR")) earlyEllis = true;
				if (peptideSequence.equals("SGGTLVLVGLGSEMTTVPLLHAAIR")) earlyEllis = true;
				if (peptideSequence.equals("SGLAACTAPLPGSR")) earlyEllis = true;
				if (peptideSequence.equals("SLSEIHCDK")) earlyEllis = true;
				if (peptideSequence.equals("VDLMAHMASKE")) earlyEllis = true;
				if (peptideSequence.equals("VYDVTQHAVGIVVNK")) earlyEllis = true;

				
				matchWriter.println("<tr>");
				
				/* counter */
				matchWriter.println("<td>" + counter + "</td>");
				
				/* acid sequence */
				matchWriter.println("<td>" + anyMatch.getString("peptideSequence") + "</td>");
				
				
				/*
				 * HACK
				 */
//				if (earlyEllis) {
//					matchWriter.println("<td>*</td>");
//				} else {
//					matchWriter.println("<td></td>");
//				}
				
				/* Gene name column */
				if (makeGeneReport) {
					int startLocus = anyMatch.getInt("start");
					int stopLocus = anyMatch.getInt("stop");
					String sequence = anyMatch.getString("sequenceName");
					if (sequence.startsWith(">")) sequence = sequence.substring(1);
					matchWriter.println("<td>");
					String interest = "";
					String geneType= "";
					for (Region region: geneRegions) {
						if (!region.getSequence().equals(sequence)) continue;
						boolean startInside = (startLocus >= region.getStart() && startLocus <= region.getStop());
						boolean stopInside = (stopLocus >= region.getStart() && stopLocus <= region.getStop());
						if (stopInside || startInside) {
							matchWriter.println("<a href=\"http://www.genecards.org/cgi-bin/carddisp.pl?gene=");
							matchWriter.println(region.getName());
							matchWriter.println("\">");
							matchWriter.println(region.getName());
							matchWriter.println("</a>");
							for (String tkName: Definitions.oncoList) {
								if (region.getName().equals(tkName)) {
									interest = "oncogene";
								}
							}
							for (String tkName: Definitions.tkList) {
								if (region.getName().equals(tkName)) {
									interest = "tyrosine-kinase";
								}
							}
							for (String tkName: Definitions.protoOncoList) {
								if (region.getName().equals(tkName)) {
									interest = "proto-oncogene";
								}
							}
							for (String tkName: Definitions.gfList) {
								if (region.getName().equals(tkName)) {
									interest = "growth factor";
								}
							}
							for (String tkName: Definitions.gtpaseActivationList) {
								if (region.getName().equals(tkName)) {
									interest = "GTPase activation";
								}
							}
							
							geneType = region.getDescription();
							break;
						}
					}
					matchWriter.println("</td>");
//					matchWriter.println("<td>" + geneType + "</td>");
					matchWriter.println("<td>" + interest + "</td>");
				}
									
				/* links */
				double scoreTotal = 0;
				int sampleCount = 0;
				int spectrumCount = 0;
				int pFiveCount = 0;
				int pSixCount = 0;
				int pFiveOrSix = 0;
				
				for (int bestMatchIndex = 0; bestMatchIndex < bestMatches.size(); bestMatchIndex++) {
					BestMatches bm = bestMatches.get(bestMatchIndex);
					Match match = holder.get(bm.getSampleName());
					
					
					if (match != null) {
						int score = (int) Math.round(match.getDouble("score"));
						
						/* if match isn't null, the peptide was present in this sample */
						sampleCount++;
						
						/* write the BED line */
						int startLocus = anyMatch.getInt("start");
						int stopLocus = anyMatch.getInt("stop");
						String sequenceName = anyMatch.getString("sequenceName");
						if (sequenceName.startsWith(">")) sequenceName = sequenceName.substring(1);
						Color color = new Color(Color.HSBtoRGB((float) bestMatchIndex / bestMatches.size(), 1.0f, 0.5f));
						String colorString = color.getRed() + "," + color.getBlue() + "," + color.getGreen();
						bedWriter.println(sequenceName + "\t" + startLocus + "\t" + stopLocus + "\t" + anyMatch.getString("peptideSequence") + "\t1000\t" + anyMatch.getString("strand") + "\t" + startLocus + "\t" + stopLocus + "\t" + colorString);
						
						/*
						 * HACK
						 */
						if (bestMatchIndex < 9) pFiveCount++;
						if (bestMatchIndex >= 9 && bestMatchIndex < 17) pSixCount++;
						if (bestMatchIndex >= 17) pFiveOrSix++;
						
						/* calculating total for all scores */
						ArrayList<Match> peptideMatches = bm.getAllPeptides().get(anyMatch.getString("peptideSequence"));
						for (Match peptideMatch: peptideMatches) {
							scoreTotal += peptideMatch.getDouble("score");
							spectrumCount++;
						}
						
							if (fullColumns) {
							/* make the spectrum link */
							matchWriter.print("<td>");
							if (makeSpectrumPage) {
	//							File spectrumFile = match.getFile("FilePath");
								File spectrumFile = new File(bm.getResultsFolder(), "spectra/" + match.getString("spectrumMD5") + ".dta");
								if (spectrumFile.exists()) {
									Spectrum spectrum = SpectrumLoader.loadSpectra(spectrumFile).get(0);
									Peptide peptide = new Peptide(match.getString("peptideSequence"));
									Match_Blank matchForReport = new Match_Blank(spectrum, peptide, score);
									ArrayList<Peppy.Match> matchesForReport = new ArrayList<Peppy.Match>();
									matchesForReport.add(matchForReport);
									
									File spectrumReportFolder = new File(parentFolder, "spectrumReportFolder");
									spectrumReportFolder.mkdir();
									File spectrumReportFile = new File(spectrumReportFolder, match.getString("spectrumMD5") + ".html");
									HTMLPageSpectrum spectrumReport = new HTMLPageSpectrum(spectrum, matchesForReport, spectrumReportFile);
									spectrumReport.makePage();
									
									/* save SVG */
									MatchSVG makeSVG = new MatchSVG(matchForReport, new File(spectrumReportFolder, match.getString("spectrumMD5") + ".svg"));
									makeSVG.saveSVG();
									
									matchWriter.print("<a href=\"spectrumReportFolder/" + spectrumReportFile.getName() + "\">");
									matchWriter.print( score );
									matchWriter.print( "</a>" );
								}
								
							} else {
								matchWriter.print( score );
							}
							matchWriter.print("</td>");
						}
						
					} else {
						if (fullColumns) {
							matchWriter.println("<td></td>");
						}
					}
				}
				
				
				
				/* score total */
				matchWriter.println("<td>" + sampleCount + "</td>");
				matchWriter.println("<td>" + spectrumCount + "</td>");
//				matchWriter.println("<td>" + holder.getScoreTotal() + "</td>");
				
				boolean isSingle = ((pFiveCount == 0) || (pSixCount == 0));
				if ( ((pFiveCount == 0) && (pSixCount == 0))) isSingle = false;
//				matchWriter.println("<td>" + pFiveCount + "</td>");
//				matchWriter.println("<td>" + pSixCount + "</td>");
//				matchWriter.println("<td>" + pFiveOrSix + "</td>");
				
//				matchWriter.println("<td>" + isSingle + "</td>");
				
				/* NIST link */
//				matchWriter.println("<td><a href=\"http://peptide.nist.gov/browser/peptide_stat.php?description=IT&organism=human&pep_seq=" + anyMatch.getString("peptideSequence") + "\">NIST</a></td>");
				
				/* sample UCSC link */
				String ucsc = UCSC.getLink(anyMatch.getInt("start"), anyMatch.getInt("stop"), anyMatch.getString("sequenceName"));
				if (showUCSC) matchWriter.println("<td><a href=\"" + ucsc + "\">UCSC</a></td>");
				
				/* unique */
				boolean unique = anyMatch.getBoolean("uniqueGlobally");
				matchWriter.println("<td>" + unique + "</td>");
				
				/* locus */
				String sequenceName = anyMatch.getString("sequenceName");
				if (sequenceName.startsWith(">")) sequenceName = sequenceName.substring(1);
				matchWriter.println("<td>" +  sequenceName + "</td>");
				matchWriter.println("<td>" + anyMatch.getInt("start") + "</td>");
				matchWriter.println("<td>" + anyMatch.getInt("stop") + "</td>");

				/* strand */
				matchWriter.println("<td>" + anyMatch.getString("strand") + "</td>");
				matchWriter.println("<td>" + note+ "</td>");
				
				matchWriter.println("</tr>");
				
				if (counter == 500) break;
				
			}
			
			U.p("in common for this report: " + counter);
			matchWriter.println("</table></ul></body></html>");

			
			/* analysis and sample legend */
			 for (int i = 0; i < bestMatches.size(); i++) {
				BestMatches bm = bestMatches.get(i);
				matchWriter.println("(" + (i + 1) + ") " + bm.getSampleName() + "; ");
			}
			
			matchWriter.flush();
			matchWriter.close();
			
			bedWriter.flush();
			bedWriter.close();

		} catch (IOException ioException) {
			ioException.printStackTrace();
		}
	}
	
	
	public void addMatchType(ResultsCategory resultsCategory) {
		this.resultsCategories.add(resultsCategory);
	}
	
	public void process() {		
		
		/* load all the matches */
		for (ResultsCategory resultsCategory: this.resultsCategories) {
			loadResults(resultsCategory);
		}
		
		populateBestPeptides();
		
		U.p("Number of spectra identified: " + bestMatches.size());
		U.p("Number of unique peptides: " + bestPeptides.size());
			
	}
	
	/**
	 * reduce best matches to best peptides
	 */
	public void populateBestPeptides() {
		bestPeptides = new Hashtable<String, Match>();
		/* reduce the best results down to the best one match for any given peptide */
		Enumeration<Match> values = bestMatches.elements();
		while (values.hasMoreElements()) {
			Match match = values.nextElement();
			String peptideSequence = match.getString("peptideSequence");
			Match bestMatch = bestPeptides.get(peptideSequence);
			if (bestMatch == null) {
				bestPeptides.put(peptideSequence, match);
				ArrayList<Match> peptideMatches = new ArrayList<Match>();
				peptideMatches.add(match);
				allPeptides.put(peptideSequence, peptideMatches);
			} else {
				
				/* add to the full list of matches for this peptide */
				ArrayList<Match> peptideMatches = allPeptides.get(peptideSequence);
				peptideMatches.add(match);
				
				/* always default to the unmodified form */
				if (match.getBoolean("isModified") == false) {
					if (bestMatch.getBoolean("isModified") == false) {
						if (match.getScore() > bestMatch.getScore() ) {
							bestPeptides.put(peptideSequence, match);
						}
					}
					
				} else {
					if (bestMatch.getBoolean("isModified")) {
						if (match.getScore() > bestMatch.getScore() ) {
							bestPeptides.put(peptideSequence, match);
						}
					} else {
						bestPeptides.put(peptideSequence, match);
					}
				}
				
			}
		}
	}
	
	
	
	private void loadResults2(ResultsCategory resultsCategory) {
		U.p("loading " + resultsCategory.getName());
		
		/* load all the matches */
		ArrayList<Match> newMatches = new ArrayList<Match>();
		
		/* this hashtable will help us determine if a match is ambiguous, 
		 * that is, the peptide (or peptide with same score)
		 * can be found elsewhere in this data set
		 */
		Hashtable<String, Match> oneMatchPerSpectrum = new Hashtable<String, Match>();
		
		for (File file: resultsCategory.getFiles()) {
			ArrayList<Match> fileMatches = Match.loadMatches(file);
			
			for (Match match: fileMatches) {
				
				/* first associate the match with this match type */
				match.set("matchType", resultsCategory);
				match.set("reportFile", file);
				
				if (match.getString("spectrumMD5") == null) U.p(file);
				Match matchToSpectrum = oneMatchPerSpectrum.get(match.getString("spectrumMD5"));
				
				
				/* all this is to determine ambiguity */
				if (matchToSpectrum == null) {
					oneMatchPerSpectrum.put(match.getString("spectrumMD5"), match);
				} else {
					if (match.getScore() > matchToSpectrum.getScore()) {
						oneMatchPerSpectrum.put(match.getString("spectrumMD5"), match);
					} else {
						if (match.getScore() == matchToSpectrum.getScore()) {
							matchToSpectrum.set("uniqueInSet",false);
							match.set("uniqueInSet",false);
						}
					}
				}
				
			}
			
			newMatches.addAll(fileMatches);
		}
		
		/* add the new results in.  This was separated out from above
		 * so that the "ambiguous" property could be accurately set
		 */
		
		for (Match match: newMatches) {
			/* Get the reigning match, if better, add */
			Match bestMatch = bestMatches.get(match.getString("spectrumMD5"));
			if (bestMatch == null) {
				bestMatches.put(match.getString("spectrumMD5"), match);
				match.set("uniqueGlobally",true);
			} else {
				if (match.getScore() > bestMatch.getScore()) {
					bestMatches.put(match.getString("spectrumMD5"), match);
					match.set("uniqueGlobally",true);
				} else {
					if (match.getScore() == bestMatch.getScore()) {
						if (match.get("amplificationScore") != null) {
							bestMatch.set("amplificationScore", match.getDouble("amplificationScore"));
						}
						bestMatch.set("uniqueGlobally", false);
						match.set("uniqueGlobally",false);
					} 
					
				}
			}
		}
	
	}

	
	
	private void loadResults(ResultsCategory resultsCategory) {
		U.p("loading " + resultsCategory.getName());
		
		for (File file: resultsCategory.getFiles()) {
			/* load one match at a time from the file */
			try {
				BufferedReader br = new BufferedReader(new FileReader(file));
				
				/* read the header lines */
				br.readLine();
				br.readLine();
				br.readLine();
				
				/* find out what property each column represents */
				String line = br.readLine();
				String [] propertyNames = line.split("\t");
				
				/* read in the first line */
				line = br.readLine();
				
				int totalSpectraAdded = 0;
				while (line != null) {
					String [] chunks = line.split("\t");
					Match match = new Match();
					for (int i = 0; i < propertyNames.length; i++) {
						Class<?> propertyType = match.getColumns().get(propertyNames[i]);
						if (propertyType == null  || propertyType.equals(String.class)) {
							match.set(propertyNames[i], chunks[i]);
						} else {
							try {
								if (propertyType.equals(File.class)) {
									match.set(propertyNames[i], new File(chunks[i]));
								}
								if (propertyType.equals(Integer.class)) {
									match.set(propertyNames[i], Integer.parseInt(chunks[i]));
								}
								if (propertyType.equals(Double.class)) {
									match.set(propertyNames[i], Double.parseDouble(chunks[i]));
								}
								if (propertyType.equals(Boolean.class)) {
									match.set(propertyNames[i], Boolean.parseBoolean(chunks[i]));
								}
							}
							catch (Exception e) {
								U.p("error with property " + propertyNames[i] + " in file " + file.getAbsolutePath());
								U.p("expected property type: " + propertyType);
								U.p("expeced property name: " + propertyNames[i]);
								U.p("column index: " + i);
								U.p("value: " + chunks[i]);
								U.p("line: " + line);
								U.p();
							}
						}
					}
					line = br.readLine();
					
					/* we now have our match */
					
					/*
					 * EXPERIMENTAL
					 */
//					if (match.getString("SequenceName").indexOf("HasPrematureStop: false") != -1) continue;
					
					/* first associate the match with this match type */
					match.set("matchType", resultsCategory);
					match.set("reportFile", file);
					
					/* Get the reigning match, if better, add */
					Match bestMatch = bestMatches.get(match.getString("spectrumMD5"));
					if (bestMatch == null) {
						bestMatches.put(match.getString("spectrumMD5"), match);
						match.set("uniqueGlobally",true);
						totalSpectraAdded++;
					} else {
						if (match.getScore() > bestMatch.getScore()) {
							bestMatches.put(match.getString("spectrumMD5"), match);
							match.set("uniqueGlobally",true);
						} else {
							if (match.getScore() == bestMatch.getScore()) {
								bestMatch.set("uniqueGlobally", false);
								match.set("uniqueGlobally",false);
							} 
							
						}
					}	
					
				}
				U.p("totalSpectraAdded: " + totalSpectraAdded);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
		}
	
	}

	/**
	 * subtracts the identified spectra from one best match set from this one
	 * @param otherBM
	 */
	public void subtractBestMatchesSpectrum(BestMatches otherBM) {
		Hashtable<String, Match> otherBestMatchesHash = otherBM.getBestMatches();
		ArrayList<String> ourBestMatchKeys = new ArrayList<String>(getBestMatches().keySet());
		Hashtable<String, Match> reducedBestMatches = new Hashtable<String, Match>();
		
		/* go through all of the keys of our best matches
		 * if the other matches don't have a key, then that means
		 * we don't lose it from the subtraction
		 */
		for (String key: ourBestMatchKeys) {
			Match otherMatch = otherBestMatchesHash.get(key);
			if (otherMatch == null) {
				reducedBestMatches.put(key, bestMatches.get(key));
			}
		}
		bestMatches = reducedBestMatches;
		populateBestPeptides();
	}
	
	/**
	 * subtracts our matches that have peptides from the set of peptides found in the other best matches
	 * @param otherBM
	 */
	public void subtractBestMatchesPeptide(BestMatches otherBM) {
		
		Hashtable<String, Match> reducedBestMatches = new Hashtable<String, Match>();
		Hashtable<String, Match> otherBestPaptides = otherBM.getBestPeptideMatches();
		ArrayList<Match> ourBestMatches = new ArrayList<Match>(getBestMatches().values());
		

		for (Match match: ourBestMatches) {
			String peptide = match.getString("peptideSequence");
			if (otherBestPaptides.get(peptide) == null) {
				reducedBestMatches.put(match.getString("spectrumMD5"), match);
			}
		}
		bestMatches = reducedBestMatches;
		populateBestPeptides();
	}
	
	
	/**
	 * makes our matches be the intersection of matches with another set
	 * @param otherBM
	 */
	public void intersectBestMatchesPeptide(BestMatches otherBM) {
		
		Hashtable<String, Match> reducedBestMatches = new Hashtable<String, Match>();
		Hashtable<String, Match> otherBestPaptides = otherBM.getBestPeptideMatches();
		ArrayList<Match> ourBestMatches = new ArrayList<Match>(getBestMatches().values());
		

		for (Match match: ourBestMatches) {
			String peptide = match.getString("peptideSequence");
			if (otherBestPaptides.get(peptide) != null) {
				reducedBestMatches.put(match.getString("spectrumMD5"), match);
			}
		}
		bestMatches = reducedBestMatches;
		populateBestPeptides();
	}


	public void saveReports() {
		ArrayList<Match> bestArray = new ArrayList<Match>(bestPeptides.values());
		Collections.sort(bestArray);
		
		File parentDirectory = new File(sampleName + " analysis");
		parentDirectory.mkdirs();
		
		/* save all matches */
		try {
			PrintWriter matchWriter = new PrintWriter(new BufferedWriter(new FileWriter(new File(parentDirectory, "all matches.txt"))));
			for (Match match: bestArray) {
				
				matchWriter.println( match.getString("peptideSequence") + "\t" + ((ResultsCategory) match.get("matchType")).getName() + "\t" + match.getFile("FilePath").getAbsolutePath());
			}
			matchWriter.flush();
			matchWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
		/* save individual match files */
		try {
			for (ResultsCategory resultsCategory: this.resultsCategories) {
				
				PrintWriter matchWriter = new PrintWriter(new BufferedWriter(new FileWriter(new File(parentDirectory, resultsCategory.getName() + " matches.html"))));
				matchWriter.println("<html>");
				matchWriter.println("<head>");
				matchWriter.println("" +
						"<style type=\"text/css\">" +
						"<!--" +
						"@import url(\"http://proteomics.me/resources/reports/style.css\");" +
						"@import url(\"http://proteomics.me/resources/reports/sortable.css\");" +
						"-->" +
						"</style>" +
						"<script type=\"text/javascript\" src=\"http://proteomics.me/resources/reports/sortable.js\"></script>"
						);
				matchWriter.println("</head>");
				matchWriter.println("<body>");
				matchWriter.println("<table class=\"sortable\" id=\"box-table-a\" >");
				matchWriter.println("<tr>");
				matchWriter.println("<th>peptide</th>");
				matchWriter.println("<th>notes</th>");
				matchWriter.println("<th>score</th>");
				if (resultsCategory.databaseType == ResultsCategory.DNA)
					matchWriter.println("<th>UCSC</th>");
				matchWriter.println("<th>NIST</th>");
				matchWriter.println("<th>sequence</th>");
				matchWriter.println("<th>start</th>");
				
				/* modification columns */
				/* isModified	modMass	modIndex	modLocCertain */
				matchWriter.println("<th>modIndex</th>");
				matchWriter.println("<th>modAcid</th>");
				matchWriter.println("<th>modMass</th>");
				matchWriter.println("</tr>");
				for (Match match: bestArray) {
					
					/*
					 * HACK
					 * for condensing the SNP_Specific reports
					 * gets rid of some uninteresting modifications
					 */
//					int modMass = (int) Math.round(match.getDouble("modMass"));
//					char modAcid = match.getString("peptideSequence").charAt(match.getInt("modIndex"));
//					if (modMass == 0) continue;
//					if (modMass == 1) continue;
//					if (modMass == 16) continue;
//					if (modMass == 144) continue;
//					if (modMass == 43) continue;
//					if (modMass == 1 && modAcid == 'N') continue;
					
					if (match.get("matchType").equals(resultsCategory)) {
						File spectrumPage = new File(match.getFile("reportFile").getParent(), "spectra/" + match.getInt("spectrumID") + ".html");
						matchWriter.println("<td><a href=\"" + spectrumPage.getAbsolutePath() + "\">" + match.getString("peptideSequence") + "</a></td>");
						matchWriter.println("<td></td>");
						matchWriter.println("<td>" + Math.round(match.getScore()) + "</td>");
						String ucsc = UCSC.getLink(match.getInt("start"), match.getInt("stop"), match.getString("SequenceName"));
						if (resultsCategory.databaseType == ResultsCategory.DNA)
							matchWriter.println("<td><a href=\"" + ucsc + "\">UCSC</a></td>");
						matchWriter.println("<td><a href=\"http://peptide.nist.gov/browser/peptide_stat.php?description=IT&organism=human&pep_seq=" + match.getString("peptideSequence") + "\">NIST</a></td>");
						if (resultsCategory.databaseType == ResultsCategory.DNA) {
							matchWriter.println("<td>" + match.getString("SequenceName") + "</td>");
						} else {
							String sequenceName = match.getString("SequenceName");
							if (sequenceName.indexOf('-') != -1) {
								sequenceName = sequenceName.substring(0, sequenceName.indexOf('-'));
							}
//							matchWriter.println("<td><a href=\"http://www.uniprot.org/uniprot/" + sequenceName + "\">" + sequenceName + "</a></td>");
							matchWriter.println("<td><a href=\"http://www.genecards.org/cgi-bin/carddisp.pl?gene=" + sequenceName + "\">" + sequenceName + "</a></td>");
							
						}
						
						
						matchWriter.println("<td>" + match.getInt("start") + "</td>");
						
						/* modification information */
//						if (match.getBoolean("isModified") && match.getBoolean("modLocCertain") ) {
							String sequence = match.getString("peptideSequence");
							int modIndex = match.getInt("modIndex");
							matchWriter.println("<td>" +  match.getInt("modIndex") + "</td>");
							matchWriter.println("<td>" + sequence.charAt(modIndex) + "</td>");
							matchWriter.println("<td>" + match.getDouble("modMass") + "</td>");
//						} else {
//							matchWriter.println("<td>" +  match.getInt("modIndex") + "</td>");
//							matchWriter.println("<td></td>");
//							matchWriter.println("<td></td>");
//						}
						
						/* end row */
						matchWriter.println("</tr>");
					}
				}
				matchWriter.println("</table></ul></body></html>");
				matchWriter.flush();
				matchWriter.close();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	
	
	public Hashtable<String, Match> getBestPeptideMatches() {
		return bestPeptides;
	}
	
	
	public Hashtable<String, Match> getBestMatches() {
		return bestMatches;
	}


	
	public String getSampleName() {
		return sampleName;
	}

	
	public ArrayList<ResultsCategory> getResultsCategories() {
		return resultsCategories;
	}


	public Hashtable<String, ArrayList<Match>> getAllPeptides() {
		return allPeptides;
	}


	public File getResultsFolder() {
		return resultsFolder;
	}


	

}
