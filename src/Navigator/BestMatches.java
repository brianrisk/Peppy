package Navigator;

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
import java.util.NoSuchElementException;

import Database.Column;
import Database.Table;
import Peppy.Properties;
import Peppy.U;
import Reports.UCSC;

public class BestMatches {
	
	/* the sample from which all of the results have been derived */
	private String sampleName;
	
	/* places on chromosomes that we want to pay attention to */
	ArrayList<Region> regionsOfInterest = new ArrayList<Region>();
	
	/* where we keep the best results */
	Hashtable<String, Match> bestMatches = new Hashtable<String, Match>();
	
	/* the results, reduced to the best example for each peptide */
	Hashtable<String, Match> bestPeptides = new Hashtable<String, Match>();
	
	/* match types */
	ArrayList<ResultsCategory> resultsCategories = new ArrayList<ResultsCategory>();
	
	public static void main(String args[]) {
//		washu();
//		washuChr8();
//		pandey();
//		mayo();
//		ucla();
		yale();
//		yaleEnzymeless();
		U.p("done");
	}
	
	public static void pandey() {
		/* Pandey */
		BestMatches pandey = new BestMatches("Pandey");
		
		/* target protein */
		ResultsCategory human = new ResultsCategory("ReferenceProtein", ResultsCategory.PROTEIN);
		human.addFile(new File("/Users/risk2/PeppyData/akhilesh-pandey/reports/Pandey/1 uncompressed - protein/report.txt"));
		pandey.addMatchType(human);
		
		/* reference genome */
		ResultsCategory genome = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
		genome.addFile(new File("/Users/risk2/PeppyData/akhilesh-pandey/reports/Pandey/2 uncompressed - genome/report.txt"));
		pandey.addMatchType(genome);
		
		/* find the best peptides */
		pandey.process();
		pandey.saveReports();
	}
	
	public static void washuChr8() {
		
		/* WHIM 16 */
		BestMatches whim16 = new BestMatches("WHIM16");
		U.p("\rloading WHIM16");
		
		/* contaminant protein */
		ResultsCategory whim16ContaminantProtein = new ResultsCategory("ContaminantProtein", ResultsCategory.PROTEIN);
		whim16ContaminantProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/chr8/WHIM16/WHIM16 nomod mouse/report.txt"));
		whim16.addMatchType(whim16ContaminantProtein);
		
		/* target protein */
		ResultsCategory whim16ReferenceProtein = new ResultsCategory("ReferenceProtein", ResultsCategory.PROTEIN);
		whim16ReferenceProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/chr8/WHIM16/WHIM16 nomod human/report.txt"));
		whim16.addMatchType(whim16ReferenceProtein);
		
		/* reference genome */
		ResultsCategory whim16ReferenceGenome = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
		whim16ReferenceGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/chr8/WHIM16/WHIM16 chr8 hg19/report.txt"));
		whim16.addMatchType(whim16ReferenceGenome);
		
		/* subject genome */
		ResultsCategory whim16SubjectGenome = new ResultsCategory("SubjectGenome", ResultsCategory.DNA);
		whim16SubjectGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/chr8/WHIM16/WHIM16 chr8 germline/report.txt"));
		whim16.addMatchType(whim16SubjectGenome);
		
		/* disease genome */
		ResultsCategory whim16DiseaseGenome = new ResultsCategory("DiseaseGenome", ResultsCategory.DNA);
		whim16DiseaseGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/chr8/WHIM16/WHIM16 chr8 xeno/report.txt"));
		whim16.addMatchType(whim16DiseaseGenome);
		
		/* find the best peptides */
		whim16.process();
		whim16.saveReports();
		
		
		
		
		/* WHIM 2 */
		BestMatches whim2 = new BestMatches("WHIM2");
		U.p("\rloading WHIM2");
		
		/* contaminant protein */
		ResultsCategory whim2ContaminantProtein = new ResultsCategory("ContaminantProtein", ResultsCategory.PROTEIN);
		whim2ContaminantProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/chr8/WHIM2/WHIM2 nomod mouse/report.txt"));
		whim2.addMatchType(whim2ContaminantProtein);
		
		/* target protein */
		ResultsCategory whim2ReferenceProtein = new ResultsCategory("ReferenceProtein", ResultsCategory.PROTEIN);
		whim2ReferenceProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/chr8/WHIM2/WHIM2 nomod human/report.txt"));
		whim2.addMatchType(whim2ReferenceProtein);
		
		/* reference genome */
		ResultsCategory whim2ReferenceGenome = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
		whim2ReferenceGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/chr8/WHIM2/WHIM2 chr8 hg19/report.txt"));
		whim2.addMatchType(whim2ReferenceGenome);
		
		/* subject genome */
		ResultsCategory whim2SubjectGenome = new ResultsCategory("SubjectGenome", ResultsCategory.DNA);
		whim2SubjectGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/chr8/WHIM2/WHIM2 chr8 germline/report.txt"));
		whim2.addMatchType(whim2SubjectGenome);
		
		/* disease genome */
		ResultsCategory whim2DiseaseGenome = new ResultsCategory("DiseaseGenome", ResultsCategory.DNA);
		whim2DiseaseGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/chr8/WHIM2/WHIM2 chr8 xeno/report.txt"));
		whim2.addMatchType(whim2DiseaseGenome);
		
		/* find the best peptides */
		whim2.process();
		whim2.saveReports();
		
		
		/* a list of our BestMatches */
		ArrayList<BestMatches> bestMatches = new ArrayList<BestMatches>();
		bestMatches.add(whim2);
		bestMatches.add(whim16);
		
		createUnifiedSamplesReport(bestMatches, "peptideSequence");
		
	}

	
public static void washuWHIM2 () {
		
//		/* WHIM 2 - 33 */
//		BestMatches whim2_33 = new BestMatches("WHIM2 - 33");
//		U.p("\rloading WHIM2 - 33");
//		
//		/* contaminant protein */
//		ResultsCategory whim2_33ContaminantProtein = new ResultsCategory("ContaminantProtein", ResultsCategory.PROTEIN);
//		whim2_33ContaminantProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/33/WHIM2/WashU WHIM2 mouse/report.txt"));
//		whim2_33.addMatchType(whim2_33ContaminantProtein);
//		
//		/* target protein */
//		ResultsCategory whim2_33ReferenceProtein = new ResultsCategory("ReferenceProtein", ResultsCategory.PROTEIN);
//		whim2_33ReferenceProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/33/WHIM2/WashU WHIM2 human/report.txt"));
//		whim2_33.addMatchType(whim2_33ReferenceProtein);
//		
//		/* reference genome */
//		ResultsCategory whim2_33ReferenceGenome = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
//		whim2_33ReferenceGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/33/WHIM2/WashU WHIM2 hg19/report.txt"));
//		whim2_33.addMatchType(whim2_33ReferenceGenome);
//		
//		/* subject genome */
//		ResultsCategory whim2_33SubjectGenome = new ResultsCategory("SubjectGenome", ResultsCategory.DNA);
//		whim2_33SubjectGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/33/WHIM2/WashU WHIM2 germline/report.txt"));
//		whim2_33.addMatchType(whim2_33SubjectGenome);
//		
//		/* disease genome */
//		ResultsCategory whim2_33DiseaseGenome = new ResultsCategory("DiseaseGenome", ResultsCategory.DNA);
//		whim2_33DiseaseGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/33/WHIM2/WashU WHIM2 xeno/report.txt"));
//		whim2_33.addMatchType(whim2_33DiseaseGenome);
//		
//		/* find the best peptides */
//		whim2_33.process();
//		whim2_33.saveReports();
		
		
		
		
		/* WHIM 2 - 41 */
		BestMatches whim2_41 = new BestMatches("WHIM2 - 41");
		U.p("\rloading WHIM2 - 41");
		
		/* contaminant protein */
		ResultsCategory whim2_41ContaminantProtein = new ResultsCategory("ContaminantProtein", ResultsCategory.PROTEIN);
		whim2_41ContaminantProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/41/WHIM2/WashU WHIM2 mouse/report.txt"));
		whim2_41.addMatchType(whim2_41ContaminantProtein);
		
		/* target protein */
		ResultsCategory whim2_41ReferenceProtein = new ResultsCategory("ReferenceProtein", ResultsCategory.PROTEIN);
		whim2_41ReferenceProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/41/WHIM2/WashU WHIM2 human/report.txt"));
		whim2_41.addMatchType(whim2_41ReferenceProtein);
		
		/* reference genome */
		ResultsCategory whim2_41ReferenceGenome = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
		whim2_41ReferenceGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/41/WHIM2/WashU WHIM2 hg19/report.txt"));
		whim2_41.addMatchType(whim2_41ReferenceGenome);
		
		/* subject genome */
		ResultsCategory whim2_41SubjectGenome = new ResultsCategory("SubjectGenome", ResultsCategory.DNA);
		whim2_41SubjectGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/41/WHIM2/WashU WHIM2 germline/report.txt"));
		whim2_41.addMatchType(whim2_41SubjectGenome);
		
		/* disease genome */
		ResultsCategory whim2_41DiseaseGenome = new ResultsCategory("DiseaseGenome", ResultsCategory.DNA);
		whim2_41DiseaseGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/41/WHIM2/WashU WHIM2 xeno/report.txt"));
		whim2_41.addMatchType(whim2_41DiseaseGenome);
		
		/* find the best peptides */
		whim2_41.process();
		whim2_41.saveReports();
		
		U.p(whim2_41.getBestMatches().get("9952a1f942ff6310afc5ff95772682fc"));
		
		
		
		
		
//		/* WHIM 2 - 43 */
//		BestMatches whim2_43 = new BestMatches("WHIM2 - 43");
//		U.p("\rloading WHIM2 - 43");
//		
//		/* contaminant protein */
//		ResultsCategory whim2_43ContaminantProtein = new ResultsCategory("ContaminantProtein", ResultsCategory.PROTEIN);
//		whim2_43ContaminantProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/43/WHIM2/WashU WHIM2 mouse/report.txt"));
//		whim2_43.addMatchType(whim2_43ContaminantProtein);
//		
//		/* target protein */
//		ResultsCategory whim2_43ReferenceProtein = new ResultsCategory("ReferenceProtein", ResultsCategory.PROTEIN);
//		whim2_43ReferenceProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/43/WHIM2/WashU WHIM2 human/report.txt"));
//		whim2_43.addMatchType(whim2_43ReferenceProtein);
//		
//		/* reference genome */
//		ResultsCategory whim2_43ReferenceGenome = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
//		whim2_43ReferenceGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/43/WHIM2/WashU WHIM2 hg19/report.txt"));
//		whim2_43.addMatchType(whim2_43ReferenceGenome);
//		
//		/* subject genome */
//		ResultsCategory whim2_43SubjectGenome = new ResultsCategory("SubjectGenome", ResultsCategory.DNA);
//		whim2_43SubjectGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/43/WHIM2/WashU WHIM2 germline/report.txt"));
//		whim2_43.addMatchType(whim2_43SubjectGenome);
//		
//		/* disease genome */
//		ResultsCategory whim2_43DiseaseGenome = new ResultsCategory("DiseaseGenome", ResultsCategory.DNA);
//		whim2_43DiseaseGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/43/WHIM2/WashU WHIM2 xeno/report.txt"));
//		whim2_43.addMatchType(whim2_43DiseaseGenome);
//		
//		/* find the best peptides */
//		whim2_43.process();
//		whim2_43.saveReports();
//		
//		
//		/* a list of our BestMatches */
		ArrayList<BestMatches> bestMatches = new ArrayList<BestMatches>();
//		bestMatches.add(whim2_33);
		bestMatches.add(whim2_41);
//		bestMatches.add(whim2_43);
//		
//		/* WHIM 16 - 33 */
//		BestMatches whim16_33 = new BestMatches("WHIM16 - 33");
//		U.p("\rloading WHIM16 - 33");
//		
//		/* contaminant protein */
//		ResultsCategory whim16_33ContaminantProtein = new ResultsCategory("ContaminantProtein", ResultsCategory.PROTEIN);
//		whim16_33ContaminantProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/33/WHIM16/WashU WHIM16 mouse/report.txt"));
//		whim16_33.addMatchType(whim16_33ContaminantProtein);
//		
//		/* target protein */
//		ResultsCategory whim16_33ReferenceProtein = new ResultsCategory("ReferenceProtein", ResultsCategory.PROTEIN);
//		whim16_33ReferenceProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/33/WHIM16/WashU WHIM16 human/report.txt"));
//		whim16_33.addMatchType(whim16_33ReferenceProtein);
//		
//		/* reference genome */
//		ResultsCategory whim16_33ReferenceGenome = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
//		whim16_33ReferenceGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/33/WHIM16/WashU WHIM16 hg19/report.txt"));
//		whim16_33.addMatchType(whim16_33ReferenceGenome);
//		
//		/* subject genome */
//		ResultsCategory whim16_33SubjectGenome = new ResultsCategory("SubjectGenome", ResultsCategory.DNA);
//		whim16_33SubjectGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/33/WHIM16/WashU WHIM16 germline/report.txt"));
//		whim16_33.addMatchType(whim16_33SubjectGenome);
//		
//		/* disease genome */
//		ResultsCategory whim16_33DiseaseGenome = new ResultsCategory("DiseaseGenome", ResultsCategory.DNA);
//		whim16_33DiseaseGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/33/WHIM16/WashU WHIM16 xeno/report.txt"));
//		whim16_33.addMatchType(whim16_33DiseaseGenome);
//		
//		/* find the best peptides */
//		whim16_33.process();
//		whim16_33.saveReports();
//		
//		
//		
//		
//		/* WHIM 16 - 41 */
//		BestMatches whim16_41 = new BestMatches("WHIM16 - 41");
//		U.p("\rloading WHIM16 - 41");
//		
//		/* contaminant protein */
//		ResultsCategory whim16_41ContaminantProtein = new ResultsCategory("ContaminantProtein", ResultsCategory.PROTEIN);
//		whim16_41ContaminantProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/41/WHIM16/WashU WHIM16 mouse/report.txt"));
//		whim16_41.addMatchType(whim16_41ContaminantProtein);
//		
//		/* target protein */
//		ResultsCategory whim16_41ReferenceProtein = new ResultsCategory("ReferenceProtein", ResultsCategory.PROTEIN);
//		whim16_41ReferenceProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/41/WHIM16/WashU WHIM16 human/report.txt"));
//		whim16_41.addMatchType(whim16_41ReferenceProtein);
//		
//		/* reference genome */
//		ResultsCategory whim16_41ReferenceGenome = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
//		whim16_41ReferenceGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/41/WHIM16/WashU WHIM16 hg19/report.txt"));
//		whim16_41.addMatchType(whim16_41ReferenceGenome);
//		
//		/* subject genome */
//		ResultsCategory whim16_41SubjectGenome = new ResultsCategory("SubjectGenome", ResultsCategory.DNA);
//		whim16_41SubjectGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/41/WHIM16/WashU WHIM16 germline/report.txt"));
//		whim16_41.addMatchType(whim16_41SubjectGenome);
//		
//		/* disease genome */
//		ResultsCategory whim16_41DiseaseGenome = new ResultsCategory("DiseaseGenome", ResultsCategory.DNA);
//		whim16_41DiseaseGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/41/WHIM16/WashU WHIM16 xeno/report.txt"));
//		whim16_41.addMatchType(whim16_41DiseaseGenome);
//		
//		/* find the best peptides */
//		whim16_41.process();
//		whim16_41.saveReports();
//		
//		
//		
//		
//		/* WHIM 16 - 43 */
//		BestMatches whim16_43 = new BestMatches("WHIM16 - 43");
//		U.p("\rloading WHIM16 - 43");
//		
//		/* contaminant protein */
//		ResultsCategory whim16_43ContaminantProtein = new ResultsCategory("ContaminantProtein", ResultsCategory.PROTEIN);
//		whim16_43ContaminantProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/43/WHIM16/WashU WHIM16 mouse/report.txt"));
//		whim16_43.addMatchType(whim16_43ContaminantProtein);
//		
//		/* target protein */
//		ResultsCategory whim16_43ReferenceProtein = new ResultsCategory("ReferenceProtein", ResultsCategory.PROTEIN);
//		whim16_43ReferenceProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/43/WHIM16/WashU WHIM16 human/report.txt"));
//		whim16_43.addMatchType(whim16_43ReferenceProtein);
//		
//		/* reference genome */
//		ResultsCategory whim16_43ReferenceGenome = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
//		whim16_43ReferenceGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/43/WHIM16/WashU WHIM16 hg19/report.txt"));
//		whim16_43.addMatchType(whim16_43ReferenceGenome);
//		
//		/* subject genome */
//		ResultsCategory whim16_43SubjectGenome = new ResultsCategory("SubjectGenome", ResultsCategory.DNA);
//		whim16_43SubjectGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/43/WHIM16/WashU WHIM16 germline/report.txt"));
//		whim16_43.addMatchType(whim16_43SubjectGenome);
//		
//		/* disease genome */
//		ResultsCategory whim16_43DiseaseGenome = new ResultsCategory("DiseaseGenome", ResultsCategory.DNA);
//		whim16_43DiseaseGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/43/WHIM16/WashU WHIM16 xeno/report.txt"));
//		whim16_43.addMatchType(whim16_43DiseaseGenome);
//		
//		/* find the best peptides */
//		whim16_43.process();
//		whim16_43.saveReports();
//		
//		
//		/* a list of our BestMatches */
//		bestMatches.add(whim16_33);
//		bestMatches.add(whim16_41);
//		bestMatches.add(whim16_43);
		
		createUnifiedSamplesReport(bestMatches, "peptideSequence");
		
		
	}
	
	public static void washuWHIM16 () {
		
		/* WHIM 16 - 33 */
		BestMatches whim16_33 = new BestMatches("WHIM16 - 33");
		U.p("\rloading WHIM16 - 33");
		
		/* contaminant protein */
		ResultsCategory whim16_33ContaminantProtein = new ResultsCategory("ContaminantProtein", ResultsCategory.PROTEIN);
		whim16_33ContaminantProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/33/WHIM16/WashU WHIM16 mouse/report.txt"));
		whim16_33.addMatchType(whim16_33ContaminantProtein);
		
		/* target protein */
		ResultsCategory whim16_33ReferenceProtein = new ResultsCategory("ReferenceProtein", ResultsCategory.PROTEIN);
		whim16_33ReferenceProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/33/WHIM16/WashU WHIM16 human/report.txt"));
		whim16_33.addMatchType(whim16_33ReferenceProtein);
		
		/* reference genome */
		ResultsCategory whim16_33ReferenceGenome = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
		whim16_33ReferenceGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/33/WHIM16/WashU WHIM16 hg19/report.txt"));
		whim16_33.addMatchType(whim16_33ReferenceGenome);
		
		/* subject genome */
		ResultsCategory whim16_33SubjectGenome = new ResultsCategory("SubjectGenome", ResultsCategory.DNA);
		whim16_33SubjectGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/33/WHIM16/WashU WHIM16 germline/report.txt"));
		whim16_33.addMatchType(whim16_33SubjectGenome);
		
		/* disease genome */
		ResultsCategory whim16_33DiseaseGenome = new ResultsCategory("DiseaseGenome", ResultsCategory.DNA);
		whim16_33DiseaseGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/33/WHIM16/WashU WHIM16 xeno/report.txt"));
		whim16_33.addMatchType(whim16_33DiseaseGenome);
		
		/* find the best peptides */
		whim16_33.process();
		whim16_33.saveReports();
		
		
		
		
		/* WHIM 16 - 41 */
		BestMatches whim16_41 = new BestMatches("WHIM16 - 41");
		U.p("\rloading WHIM16 - 41");
		
		/* contaminant protein */
		ResultsCategory whim16_41ContaminantProtein = new ResultsCategory("ContaminantProtein", ResultsCategory.PROTEIN);
		whim16_41ContaminantProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/41/WHIM16/WashU WHIM16 mouse/report.txt"));
		whim16_41.addMatchType(whim16_41ContaminantProtein);
		
		/* target protein */
		ResultsCategory whim16_41ReferenceProtein = new ResultsCategory("ReferenceProtein", ResultsCategory.PROTEIN);
		whim16_41ReferenceProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/41/WHIM16/WashU WHIM16 human/report.txt"));
		whim16_41.addMatchType(whim16_41ReferenceProtein);
		
		/* reference genome */
		ResultsCategory whim16_41ReferenceGenome = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
		whim16_41ReferenceGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/41/WHIM16/WashU WHIM16 hg19/report.txt"));
		whim16_41.addMatchType(whim16_41ReferenceGenome);
		
		/* subject genome */
		ResultsCategory whim16_41SubjectGenome = new ResultsCategory("SubjectGenome", ResultsCategory.DNA);
		whim16_41SubjectGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/41/WHIM16/WashU WHIM16 germline/report.txt"));
		whim16_41.addMatchType(whim16_41SubjectGenome);
		
		/* disease genome */
		ResultsCategory whim16_41DiseaseGenome = new ResultsCategory("DiseaseGenome", ResultsCategory.DNA);
		whim16_41DiseaseGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/41/WHIM16/WashU WHIM16 xeno/report.txt"));
		whim16_41.addMatchType(whim16_41DiseaseGenome);
		
		/* find the best peptides */
		whim16_41.process();
		whim16_41.saveReports();
		
		
		
		
		/* WHIM 16 - 43 */
		BestMatches whim16_43 = new BestMatches("WHIM16 - 43");
		U.p("\rloading WHIM16 - 43");
		
		/* contaminant protein */
		ResultsCategory whim16_43ContaminantProtein = new ResultsCategory("ContaminantProtein", ResultsCategory.PROTEIN);
		whim16_43ContaminantProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/43/WHIM16/WashU WHIM16 mouse/report.txt"));
		whim16_43.addMatchType(whim16_43ContaminantProtein);
		
		/* target protein */
		ResultsCategory whim16_43ReferenceProtein = new ResultsCategory("ReferenceProtein", ResultsCategory.PROTEIN);
		whim16_43ReferenceProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/43/WHIM16/WashU WHIM16 human/report.txt"));
		whim16_43.addMatchType(whim16_43ReferenceProtein);
		
		/* reference genome */
		ResultsCategory whim16_43ReferenceGenome = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
		whim16_43ReferenceGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/43/WHIM16/WashU WHIM16 hg19/report.txt"));
		whim16_43.addMatchType(whim16_43ReferenceGenome);
		
		/* subject genome */
		ResultsCategory whim16_43SubjectGenome = new ResultsCategory("SubjectGenome", ResultsCategory.DNA);
		whim16_43SubjectGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/43/WHIM16/WashU WHIM16 germline/report.txt"));
		whim16_43.addMatchType(whim16_43SubjectGenome);
		
		/* disease genome */
		ResultsCategory whim16_43DiseaseGenome = new ResultsCategory("DiseaseGenome", ResultsCategory.DNA);
		whim16_43DiseaseGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/43/WHIM16/WashU WHIM16 xeno/report.txt"));
		whim16_43.addMatchType(whim16_43DiseaseGenome);
		
		/* find the best peptides */
		whim16_43.process();
		whim16_43.saveReports();
		
		
		/* a list of our BestMatches */
		ArrayList<BestMatches> bestMatches = new ArrayList<BestMatches>();
		bestMatches.add(whim16_33);
		bestMatches.add(whim16_41);
		bestMatches.add(whim16_43);
		
		createUnifiedSamplesReport(bestMatches, "peptideSequence");
		
		
	}
	
	public static void washu41() {
//		/* WHIM 2 */
//		BestMatches whim2 = new BestMatches("WHIM2");
//		U.p("\rloading WHIM2");
//		
//		/* contaminant protein */
//		ResultsCategory whim2ContaminantProtein = new ResultsCategory("ContaminantProtein", ResultsCategory.PROTEIN);
//		whim2ContaminantProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/WHIM2/WashU WHIM2 mouse/report.txt"));
//		whim2.addMatchType(whim2ContaminantProtein);
//		
//		/* target protein */
//		ResultsCategory whim2ReferenceProtein = new ResultsCategory("ReferenceProtein", ResultsCategory.PROTEIN);
//		whim2ReferenceProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/WHIM2/WashU WHIM2 human/report.txt"));
//		whim2.addMatchType(whim2ReferenceProtein);
//		
//		/* reference genome */
//		ResultsCategory whim2ReferenceGenome = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
//		whim2ReferenceGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/WHIM2/WashU WHIM2 hg19/report.txt"));
//		whim2.addMatchType(whim2ReferenceGenome);
//		
//		/* subject genome */
//		ResultsCategory whim2SubjectGenome = new ResultsCategory("SubjectGenome", ResultsCategory.DNA);
//		whim2SubjectGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/WHIM2/WashU WHIM2 germline/report.txt"));
//		whim2.addMatchType(whim2SubjectGenome);
//		
//		/* disease genome */
//		ResultsCategory whim2DiseaseGenome = new ResultsCategory("DiseaseGenome", ResultsCategory.DNA);
//		whim2DiseaseGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/WHIM2/WashU WHIM2 xeno/report.txt"));
//		whim2.addMatchType(whim2DiseaseGenome);
//		
//		
//		/* find the best peptides */
//		whim2.process();
//		whim2.saveReports();
		
		
		/* this list should be ordered by hierarchy of importance; least important first */
		
		/* WHIM 16 */
		BestMatches whim16 = new BestMatches("WHIM16");
		U.p("\rloading WHIM16");
		
		/* contaminant protein */
		ResultsCategory whim16ContaminantProtein = new ResultsCategory("ContaminantProtein", ResultsCategory.PROTEIN);
		whim16ContaminantProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/41/WHIM16/WashU WHIM16 mouse/report.txt"));
		whim16.addMatchType(whim16ContaminantProtein);
		
		/* target protein */
		ResultsCategory whim16ReferenceProtein = new ResultsCategory("ReferenceProtein", ResultsCategory.PROTEIN);
		whim16ReferenceProtein.addFile(new File("/Users/risk2/PeppyData/WashU/reports/41/WHIM16/WashU WHIM16 human/report.txt"));
		whim16.addMatchType(whim16ReferenceProtein);
		
		/* reference genome */
		ResultsCategory whim16ReferenceGenome = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
		whim16ReferenceGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/41/WHIM16/WashU WHIM16 hg19/report.txt"));
		whim16.addMatchType(whim16ReferenceGenome);
		
		/* subject genome */
		ResultsCategory whim16SubjectGenome = new ResultsCategory("SubjectGenome", ResultsCategory.DNA);
		whim16SubjectGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/41/WHIM16/WashU WHIM16 germline/report.txt"));
		whim16.addMatchType(whim16SubjectGenome);
		
		/* disease genome */
		ResultsCategory whim16DiseaseGenome = new ResultsCategory("DiseaseGenome", ResultsCategory.DNA);
		whim16DiseaseGenome.addFile(new File("/Users/risk2/PeppyData/WashU/reports/41/WHIM16/WashU WHIM16 xeno/report.txt"));
		whim16.addMatchType(whim16DiseaseGenome);
		
		/* find the best peptides */
		whim16.process();
		whim16.saveReports();
		
		
		/* a list of our BestMatches */
		ArrayList<BestMatches> bestMatches = new ArrayList<BestMatches>();
//		bestMatches.add(whim2);
		bestMatches.add(whim16);
		
		createUnifiedSamplesReport(bestMatches, "peptideSequence");
		
		
	}

	

	
	public static void washu() {
		/* WHIM 2 */
		BestMatches whim2 = new BestMatches("WHIM2");
		U.p("\rloading WHIM2");
		
		/* contaminant protein */
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
	
	public static void yale() {
		Properties.isYale = true;
		
		
		/* Yale 2011-06 */
		BestMatches mouse1 = new BestMatches("2011-06");
		
		/* contaminant protein */
		ResultsCategory mouse1Contaminant = new ResultsCategory("ContaminantProtein", ResultsCategory.PROTEIN);
		mouse1Contaminant.addFile(new File("/Users/risk2/Sites/research/karen-anderson/yale-mouse 2011-06/yale-mouse 2011-06_report.txt"));
		mouse1.addMatchType(mouse1Contaminant);

		/* reference genome */
		ResultsCategory mouse1Genome = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
		mouse1Genome.addFile(new File("/Users/risk2/Sites/research/karen-anderson/yale-genome 2011-06/yale-genome 2011-06_report.txt"));
		mouse1.addMatchType(mouse1Genome);
		
		
		/* find the best peptides */
		mouse1.process();
		mouse1.saveReports();
		
		
		
		/* Yale 2011-10 */
		BestMatches mouse2 = new BestMatches("2011-10");
		
		/* contaminant protein */
		ResultsCategory mouse2Contaminant = new ResultsCategory("ContaminantProtein", ResultsCategory.PROTEIN);
		mouse2Contaminant.addFile(new File("/Users/risk2/Sites/research/karen-anderson/yale-mouse 2011-10/yale-mouse 2011-10_report.txt"));
		mouse2.addMatchType(mouse2Contaminant);

		/* reference genome */
		ResultsCategory mouse2Reference = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
		mouse2Reference.addFile(new File("/Users/risk2/Sites/research/karen-anderson/yale-genome 2011-10/yale-genome 2011-10_report.txt"));
		mouse2.addMatchType(mouse2Reference);
		
		
		/* find the best peptides */
		mouse2.process();
		mouse2.saveReports();
		
		
		
		/* Yale 2011-11 */
		BestMatches mouse3 = new BestMatches("2011-11");
		
		/* contaminant protein */
		ResultsCategory mouse3Contaminant = new ResultsCategory("ContaminantProtein", ResultsCategory.PROTEIN);
		mouse3Contaminant.addFile(new File("/Users/risk2/Sites/research/karen-anderson/yale-mouse 2011-11/yale-mouse 2011-11_report.txt"));
		mouse3.addMatchType(mouse3Contaminant);

		/* reference genome */
		ResultsCategory mouse3Reference = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
		mouse3Reference.addFile(new File("/Users/risk2/Sites/research/karen-anderson/yale-genome 2011-11/yale-genome 2011-11_report.txt"));
		mouse3.addMatchType(mouse3Reference);
		
		
		/* find the best peptides */
		mouse3.process();
		mouse3.saveReports();
		
		
		
		
		
		/* Yale 2012-04 tryptic*/
		BestMatches mouse4A = new BestMatches("04 tryptic");
		
		/* contaminant protein */
		ResultsCategory mouse4AContaminant = new ResultsCategory("ContaminantProtein", ResultsCategory.PROTEIN);
		mouse4AContaminant.addFile(new File("/Users/risk2/Sites/research/karen-anderson/2012-06-02/2012-04 trypsin trypsin/1 trypsin - MOUSE.fasta/report.txt"));
		mouse4A.addMatchType(mouse4AContaminant);

		/* reference genome */
		ResultsCategory mouse4AReference = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
		mouse4AReference.addFile(new File("/Users/risk2/Sites/research/karen-anderson/2012-06-02/2012-04 trypsin trypsin/2 trypsin - mouse/report.txt"));
		mouse4A.addMatchType(mouse4AReference);
		
		
		/* find the best peptides */
		mouse4A.process();
		mouse4A.saveReports();
		
		/* Yale 2012-04 chymotryptic*/
		BestMatches mouse4B = new BestMatches("04 chymo");
		
		/* contaminant protein */
		ResultsCategory mouse4BContaminant = new ResultsCategory("ContaminantProtein", ResultsCategory.PROTEIN);
		mouse4BContaminant.addFile(new File("/Users/risk2/Sites/research/karen-anderson/2012-06-02/2012-04 chymo chymo/1 chymotrypsin - MOUSE.fasta/report.txt"));
		mouse4B.addMatchType(mouse4BContaminant);

		/* reference genome */
		ResultsCategory mouse4BReference = new ResultsCategory("ReferenceGenome", ResultsCategory.DNA);
		mouse4BReference.addFile(new File("/Users/risk2/Sites/research/karen-anderson/2012-06-02/2012-04 chymo chymo/2 chymotrypsin - mouse/report.txt"));
		mouse4B.addMatchType(mouse4BReference);
		
		
		/* find the best peptides */
		mouse4B.process();
		mouse4B.saveReports();
		
		

		
		
		/* a list of our BestMatches */
		ArrayList<BestMatches> bestMatches = new ArrayList<BestMatches>();
		bestMatches.add(mouse1);
		bestMatches.add(mouse2);
		bestMatches.add(mouse3);
		bestMatches.add(mouse4A);
		bestMatches.add(mouse4B);
		
		
		createUnifiedSamplesReport(bestMatches, "peptideSequence");
		
	}
	
	
	public static void yaleEnzymeless() {
		Properties.isYale = true;
		
		
		/* Yale 2011-06 */
		BestMatches enzymeless = new BestMatches("2012-04 enzymeless");
		
		
		
		/* tryp trypsin*/
		ResultsCategory trypsinTrypsin = new ResultsCategory("trypsinTrypsin", ResultsCategory.PROTEIN);
		trypsinTrypsin.addFile(new File("/Users/risk2/Sites/research/karen-anderson/2012-06-02/2012-04 trypsin trypsin/1 trypsin - MOUSE.fasta/report.txt"));
		enzymeless.addMatchType(trypsinTrypsin);

		/* chymo chymo */
		ResultsCategory chymoChymo = new ResultsCategory("chymoChymo", ResultsCategory.PROTEIN);
		chymoChymo.addFile(new File("/Users/risk2/Sites/research/karen-anderson/2012-06-02/2012-04 chymo chymo/2 chymotrypsin - mouse/report.txt"));
		enzymeless.addMatchType(chymoChymo);
		
		
		/* enzymeless */
		ResultsCategory noEnzyeReulsts = new ResultsCategory("enzymeless", ResultsCategory.PROTEIN);
		noEnzyeReulsts.addFile(new File("/Users/risk2/Sites/research/karen-anderson/2012-06-02/2012-04 a enzymeless/1 2012-04 - MOUSE.fasta/report.txt"));
		enzymeless.addMatchType(noEnzyeReulsts);
		
		
		
		
		
		/* find the best peptides */
		enzymeless.process();
		enzymeless.saveReports();
	}
	
	
	public BestMatches(String sampleName) {
		this.sampleName = sampleName;
	}
	
	public static void createUnifiedSamplesReport(ArrayList<BestMatches> bestMatches, String keyName) {
		/* label the matches for their respective samples */
		Enumeration<Match> e;
		for (BestMatches bm: bestMatches) {
			e = bm.getBestMatches().elements();
			while (e.hasMoreElements()) e.nextElement().set("sampleName", bm.getSampleName());
		}

		
		/* create tables of results for merging */
		Column key = new Column(keyName, String.class);
		ArrayList<Table> tables = new ArrayList<Table>(bestMatches.size());
		for (BestMatches bm: bestMatches) {
			Table table = new Table(bm.getSampleName(), key);
			e = bm.getBestMatches().elements();
			while (e.hasMoreElements()) table.add(e.nextElement());
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
		
		
		
		
		/* write merged table */
		try {
			PrintWriter matchWriter = new PrintWriter(new BufferedWriter(new FileWriter(new File("unified samples report.html"))));
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
			matchWriter.println("<th>#</th>");
			matchWriter.println("<th>peptide</th>");
			for (BestMatches bm: bestMatches) {
				matchWriter.println("<th>" + bm.getSampleName() + "</th>");
			}
			matchWriter.println("<th>modified</th>");
			matchWriter.println("<th>UCSC</th>");
			matchWriter.println("<th>NIST</th>");
			matchWriter.println("</tr>");
			
			int counter = 0;
			for (MergeMatchHolder holder: pmmhs) {
				Match anyMatch;
				anyMatch = holder.get();
				if (((ResultsCategory) anyMatch.get("matchType")).databaseType == ResultsCategory.DNA) {
//				if (holder.get("ucla loo deprived") != null ) {
//					if(holder.get("ucla loo deprived").getBoolean("isModified")) {
					
					counter++;
					
					/* counter */
					matchWriter.println("<td>" + counter + "</td>");
					
					/* acid sequence */
					matchWriter.println("<td>" + anyMatch.getString("peptideSequence") + "</td>");
					
					/* links */
					for (BestMatches bm: bestMatches) {
						Match match = holder.get(bm.getSampleName());
						if (match != null) {
							File spectrumPage = new File(match.getFile("reportFile").getParent(), "spectra/" + match.getInt("spectrumID") + ".html");
							int score = (int) Math.round(match.getDouble("score"));
							matchWriter.print("<td>");
							matchWriter.print("<a href=\"" + spectrumPage.getAbsolutePath());
							matchWriter.print("\">" + score + "</a>");
							if (match.getBoolean("isModified")) matchWriter.print(" (M)");
							matchWriter.print("</td>");
						} else {
							matchWriter.println("<td></td>");
						}
					}
					
					

					
					
					/* if modified*/
					matchWriter.println("<td>" + holder.isModified() + "</td>");
					
					/* sample UCSC link */
					String ucsc = UCSC.getLink(anyMatch.getInt("start"), anyMatch.getInt("stop"), anyMatch.getString("sequenceName"));
					matchWriter.println("<td><a href=\"" + ucsc + "\">UCSC</a></td>");
					
					/* NIST link */
					matchWriter.println("<td><a href=\"http://peptide.nist.gov/browser/peptide_stat.php?description=IT&organism=human&pep_seq=" + anyMatch.getString("peptideSequence") + "\">NIST</a></td>");
					
					matchWriter.println("</tr>");
					
				}
			}
			
			U.p("in common for this report: " + counter);
			matchWriter.println("</table></ul></body></html>");
			matchWriter.flush();
			matchWriter.close();

		} catch (IOException ioException) {
			ioException.printStackTrace();
		}
	}
	
	
	public void addMatchType(ResultsCategory resultsCategory) {
		this.resultsCategories.add(resultsCategory);
	}
	
	public void process() {
		/* load the regions of interest */
//		loadRegionsOfInterest(new File("/Users/risk2/Documents/WashU/aCGH_whim2_hg19_all.txt"));

		
		
		/* load all the matches */
		for (ResultsCategory resultsCategory: this.resultsCategories) {
			loadResults(resultsCategory);
		}
				
		
		/* reduce the best results down to the best one match for any given peptide */
		Enumeration<Match> values = bestMatches.elements();
		while (values.hasMoreElements()) {
			Match match = values.nextElement();
			String peptideSequence = match.getString("peptideSequence");
			Match bestMatch = bestPeptides.get(peptideSequence);
			if (bestMatch == null) {
				bestPeptides.put(peptideSequence, match);
			} else {
				
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
		
		U.p("Number of spectra identified: " + bestMatches.size());
		U.p("Number of unique peptides: " + bestPeptides.size());
			
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
				matchWriter.println("<th>score</th>");
				matchWriter.println("<th>modified</th>");
				if (resultsCategory.databaseType == ResultsCategory.DNA)
					matchWriter.println("<th>UCSC</th>");
				matchWriter.println("<th>NIST</th>");
				matchWriter.println("<th>sequence</th>");
				matchWriter.println("<th>start</th>");
				matchWriter.println("</tr>");
				for (Match match: bestArray) {
					if (match.get("matchType").equals(resultsCategory)) {
						File spectrumPage = new File(match.getFile("reportFile").getParent(), "spectra/" + match.getInt("spectrumID") + ".html");
						matchWriter.println("<td><a href=\"" + spectrumPage.getAbsolutePath() + "\">" + match.getString("peptideSequence") + "</a></td>");
						matchWriter.println("<td>" + Math.round(match.getScore()) + "</td>");
						matchWriter.println("<td>" + match.getBoolean("isModified")+ "</td>");
						String ucsc = UCSC.getLink(match.getInt("start"), match.getInt("stop"), match.getString("SequenceName"));
						if (resultsCategory.databaseType == ResultsCategory.DNA)
							matchWriter.println("<td><a href=\"" + ucsc + "\">UCSC</a></td>");
						matchWriter.println("<td><a href=\"http://peptide.nist.gov/browser/peptide_stat.php?description=IT&organism=human&pep_seq=" + match.getString("peptideSequence") + "\">NIST</a></td>");
						if (resultsCategory.databaseType == ResultsCategory.DNA) {
							matchWriter.println("<td>" + match.getString("SequenceName") + "</td>");
						} else {
							String sequenceName = match.getString("SequenceName");
							matchWriter.println("<td><a href=\"http://www.uniprot.org/uniprot/" + sequenceName + "\">" + sequenceName + "</a></td>");
						}
						
						
						matchWriter.println("<td>" + match.getInt("start") + "</td>");
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
	
	
	
	private void loadRegionsOfInterest(File regionsFile) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(regionsFile));
			
			String line = br.readLine();
			
			while (line != null) {
				String [] chunks = line.split("\t");
				regionsOfInterest.add(
						new Region(
							chunks[0],
							Integer.parseInt(chunks[1]),
							Integer.parseInt(chunks[2]),
							Double.parseDouble(chunks[3])
						)
					);
				/* don't forget to load in that next line! */
				line = br.readLine();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}



	
	private void loadResults(ResultsCategory resultsCategory) {
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
				
				if (match.getString("peptideSequence").equals("GQQALIASSGLTVEVDAPK")) U.p(match);
				
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
				
				/* determine if match is in a region of interest and set that property */
				if (resultsCategory.getDatabaseType() == resultsCategory.DNA) {
					for (Region region: regionsOfInterest) {
						if (region.getStart() <= match.getInt("start")) {
							if (region.getStop() >= match.getInt("stop")) {
								if (region.getSequenceName().equalsIgnoreCase(match.getString("sequenceName"))) {
									match.set("amplificationScore", region.getScore());
								}
							}
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
			} else {
				if (match.getScore() > bestMatch.getScore()) {
					bestMatches.put(match.getString("spectrumMD5"), match);
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
	
	
	public String getSampleName() {
		return sampleName;
	}

	
	public ArrayList<ResultsCategory> getResultsCategories() {
		return resultsCategories;
	}
	
	

}
