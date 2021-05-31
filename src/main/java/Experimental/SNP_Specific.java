package Experimental;

import Navigator.BestMatches;
import Navigator.Match;
import Navigator.ResultsCategory;
import Peppy.U;

import java.io.File;
import java.util.ArrayList;

public class SNP_Specific {
	
	public static void main(String [] args) {
//		File modResultsFile = new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM2-Ellis043/5 WHIM2 - varimod/report.txt");
//		File snpResultsFile = new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM16-Ellis043/2 WHIM16 - subject/report.txt");
//		File modResultsFile = new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM16-Ellis043/6 WHIM16 - varimod/report.txt");
//		File snpResultsFile = new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM2-Ellis043/2 WHIM2 - subject/report.txt");
//		snpRemovesModSite(modResultsFile, snpResultsFile);
		
//		snpSpecificPTM();
		
		pccSNPSpecific();
		
		U.p("done");
		
	}
	
	public static void snpSpecificPTM() {
//		BestMatches whim2 = new BestMatches("WHIM2");
//		ResultsCategory results = new ResultsCategory("WHIM 2 SNP-specific matches", ResultsCategory.PROTEIN);
//		results.addFile(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM2-Ellis033/2 WHIM2 - subject/report.txt"));
//		results.addFile(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM2-Ellis041/2 WHIM2 - subject/report.txt"));
//		results.addFile(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM2-Ellis043/2 WHIM2 - subject/report.txt"));
//		results.addFile(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM2-Ellis043/2 WHIM2 - subject/report.txt"));
//		whim2.addMatchType(results);
//		whim2.process();
//		
//		BestMatches whim2Varimod = new BestMatches("WHIM2 varimod");
//		ResultsCategory resultsVariMod = new ResultsCategory("WHIM 2 SNP-specific matches", ResultsCategory.PROTEIN);
//		resultsVariMod.addFile(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM2-Ellis033/5 WHIM2 - varimod/report.txt"));
//		resultsVariMod.addFile(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM2-Ellis041/6 WHIM2 - varimod/report.txt"));
//		resultsVariMod.addFile(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM2-Ellis043/5 WHIM2 - varimod/report.txt"));
//		resultsVariMod.addFile(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM2-Ellis043/5 WHIM2 - varimod/report.txt"));
//		whim2Varimod.addMatchType(resultsVariMod);
//		whim2Varimod.process();
//		
//		whim2Varimod.intersectBestMatchesPeptide(whim2);
//		whim2Varimod.saveReports();
//		
		
//		BestMatches whim16 = new BestMatches("WHIM16");
//		ResultsCategory resultsWhim16 = new ResultsCategory("WHIM 16 SNP-specific matches", ResultsCategory.PROTEIN);
//		resultsWhim16.addFile(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM16-Ellis033/2 WHIM16 - subject/report.txt"));
//		resultsWhim16.addFile(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM16-Ellis041/2 WHIM16 - subject/report.txt"));
//		resultsWhim16.addFile(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM16-Ellis043/2 WHIM16 - subject/report.txt"));
//		resultsWhim16.addFile(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM16-Ellis043/2 WHIM16 - subject/report.txt"));
//		whim16.addMatchType(resultsWhim16);
//		whim16.process();
//		
//		BestMatches whim16Varimod = new BestMatches("WHIM16 varimod");
//		ResultsCategory resultsWhim16VariMod = new ResultsCategory("WHIM 16 SNP-specific matches", ResultsCategory.PROTEIN);
//		resultsWhim16VariMod.addFile(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM16-Ellis033/6 WHIM16 - varimod/report.txt"));
//		resultsWhim16VariMod.addFile(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM16-Ellis041/6 WHIM16 - varimod/report.txt"));
//		resultsWhim16VariMod.addFile(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM16-Ellis043/5 WHIM16 - varimod/report.txt"));
//		resultsWhim16VariMod.addFile(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM16-Ellis043/6 WHIM16 - varimod/report.txt"));
//		whim16Varimod.addMatchType(resultsWhim16VariMod);
//		whim16Varimod.process();
//		
//		whim16Varimod.intersectBestMatchesPeptide(whim16);
//		whim16Varimod.saveReports();
		
		
//		BestMatches whim16 = new BestMatches("WHIM16");
//		ResultsCategory resultsWhim16 = new ResultsCategory("WHIM 16 SNP-specific matches", ResultsCategory.PROTEIN);
//		resultsWhim16.addFile(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/UNC-WHIM16-CompRef/2 WHIM16 - germline proteome/report.txt"));
//
//		whim16.addMatchType(resultsWhim16);
//		whim16.process();
//		
//		BestMatches whim16Varimod = new BestMatches("WHIM16 varimod");
//		ResultsCategory resultsWhim16VariMod = new ResultsCategory("WHIM 16 SNP-specific matches", ResultsCategory.PROTEIN);
//		resultsWhim16VariMod.addFile(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/UNC-WHIM16-CompRef/6 WHIM16 - varimod/report.txt"));
//
//		whim16Varimod.addMatchType(resultsWhim16VariMod);
//		whim16Varimod.process();
//		
//		whim16Varimod.intersectBestMatchesPeptide(whim16);
//		whim16Varimod.saveReports();
		
		
		
//		BestMatches whim2 = new BestMatches("WHIM2");
//		ResultsCategory resultsWhim2 = new ResultsCategory("WHIM 2 SNP-specific matches", ResultsCategory.PROTEIN);
//		resultsWhim2.addFile(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/UNC-WHIM2-CompRef/3 WHIM2 - xeno proteome/report.txt"));
//
//		whim2.addMatchType(resultsWhim2);
//		whim2.process();
//		
//		BestMatches whim2Varimod = new BestMatches("WHIM2 varimod");
//		ResultsCategory resultsWhim2VariMod = new ResultsCategory("WHIM 2 SNP-specific matches", ResultsCategory.PROTEIN);
//		resultsWhim2VariMod.addFile(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/UNC-WHIM2-CompRef/6 WHIM2 - varimod/report.txt"));
//
//		whim2Varimod.addMatchType(resultsWhim2VariMod);
//		whim2Varimod.process();
//		
//		whim2Varimod.intersectBestMatchesPeptide(whim2);
//		whim2Varimod.saveReports();
		
		BestMatches whim16 = new BestMatches("WHIM16");
		ResultsCategory resultsWhim16 = new ResultsCategory("WHIM 16 SNP-specific matches", ResultsCategory.PROTEIN);
		resultsWhim16.addFile(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/UNC-WHIM16-CompRef/3 WHIM16 - xeno proteome/report.txt"));

		whim16.addMatchType(resultsWhim16);
		whim16.process();
		
		BestMatches whim16Varimod = new BestMatches("WHIM16 varimod");
		ResultsCategory resultsWhim16VariMod = new ResultsCategory("WHIM 16 SNP-specific matches", ResultsCategory.PROTEIN);
		resultsWhim16VariMod.addFile(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/UNC-WHIM16-CompRef/6 WHIM16 - varimod/report.txt"));

		whim16Varimod.addMatchType(resultsWhim16VariMod);
		whim16Varimod.process();
		
		whim16Varimod.intersectBestMatchesPeptide(whim16);
		whim16Varimod.saveReports();
		
		
		
	
		

		
		//GVVDSEEIPLNLSR
		
	}
	
	public static void pccSNPSpecific() {

		ArrayList<File> reportFolders = new ArrayList<File>();
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM2-Ellis041/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM2-Ellis043/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM2-Ellis043/"));
		reportFolders.add(new File("reports/CompRef_Proteome_VU_B1_P5"));
		reportFolders.add(new File("reports/CompRef_Proteome_VU_B2_P5"));
		reportFolders.add(new File("reports/CompRef_Proteome_VU_B3_P5"));
		reportFolders.add(new File("reports/CompRef_Proteome_JHUC_P5AB"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM2-CompRef/"));
		
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM16-Ellis041/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM16-Ellis043/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM16-Ellis043/"));
		reportFolders.add(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM16-CompRef/"));
		reportFolders.add(new File("reports/CompRef_Proteome_VU_A2_P6"));
		reportFolders.add(new File("reports/CompRef_Proteome_VU_A3_P6"));		
		reportFolders.add(new File("reports/CompRef_Proteome_JHUC_P6ST"));
		
		
		reportFolders.add(new File("reports/CompRef_Proteome_PNNL"));
		reportFolders.add(new File("reports/CompRef_Proteome_JHUC_iTRAQ"));
		reportFolders.add(new File("reports/CompRef_Proteome_BI"));
		
		/* a list of folders to ignore from our results */
		ArrayList<String> direcotryTitlesToIgnoreForNoMods = new ArrayList<String>();
		direcotryTitlesToIgnoreForNoMods.add("personal");
//		direcotryTitlesToIgnoreForNoMods.add("subject");
		direcotryTitlesToIgnoreForNoMods.add("varimod");
		direcotryTitlesToIgnoreForNoMods.add("HG19");
		direcotryTitlesToIgnoreForNoMods.add("mouse");
//		direcotryTitlesToIgnoreForNoMods.add("germline");
		direcotryTitlesToIgnoreForNoMods.add("xeno");
		direcotryTitlesToIgnoreForNoMods.add("gencode");
		
		BestMatches matches = new BestMatches(new File("reports/broad"), -1, direcotryTitlesToIgnoreForNoMods);
		
		/* a list of folders to ignore from our results */
		ArrayList<String> direcotryTitlesToIgnoreForMods = new ArrayList<String>();
		direcotryTitlesToIgnoreForMods.add("personal");
		direcotryTitlesToIgnoreForMods.add("subject");
//		direcotryTitlesToIgnoreForMods.add("varimod");
		direcotryTitlesToIgnoreForMods.add("HG19");
		direcotryTitlesToIgnoreForMods.add("mouse");
		direcotryTitlesToIgnoreForMods.add("germline");
		direcotryTitlesToIgnoreForMods.add("xeno");
		direcotryTitlesToIgnoreForMods.add("gencode");
		
		
		for (File reportFolder: reportFolders) {
			BestMatches modMatches = new BestMatches(reportFolder, -1, direcotryTitlesToIgnoreForMods);
			modMatches.intersectBestMatchesPeptide(matches);
			modMatches.saveReports();
		}
		
		
		
	}
	
	/**
	 * This compares two different spectral sets against each other.
	 * It finds that if a nsSNP is in one, and in the other there is a mod
	 * where the nsSNP is and (of course) the AA is different, then that implies
	 * the nsSNP has removed that mod site.
	 */
	public static void snpRemovesModSite(File modResultsFile, File snpResultsFile) {
		ArrayList<Match> modResults = Match.loadMatches(modResultsFile);
		modResults.addAll(Match.loadMatches(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM2-Ellis033/5 WHIM2 - varimod/report.txt")));
		modResults.addAll(Match.loadMatches(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM2-Ellis041/6 WHIM2 - varimod/report.txt")));
		modResults.addAll(Match.loadMatches(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM2-Ellis043/5 WHIM2 - varimod/report.txt")));
		
		ArrayList<Match> snpResults = Match.loadMatches(snpResultsFile);
		snpResults.addAll(Match.loadMatches(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM16-Ellis033/2 WHIM16 - subject/report.txt")));
		snpResults.addAll(Match.loadMatches(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM16-Ellis041/2 WHIM16 - subject/report.txt")));
		snpResults.addAll(Match.loadMatches(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM16-Ellis043/2 WHIM16 - subject/report.txt")));
		
//		ArrayList<Match> modResultsThat
		int differenceLoc;
		for (Match modMatch: modResults) {
			int roundedMass = (int) Math.round(modMatch.getDouble("modMass"));
			
			/*
			 * skipping non-interesting mods
			 */
			if (roundedMass == 1) continue; 
			if (roundedMass == 2) continue; 
			if (roundedMass == 43) continue; 
			if (roundedMass == 16) continue;
			
			for (Match snpMatch: snpResults) {
				differenceLoc = oneDifference(modMatch,snpMatch);
				if (differenceLoc != -1) {
					if (modMatch.getInt("modIndex") == differenceLoc) {
						U.p(roundedMass + "\t" + modMatch.getString("peptideSequence") + "\t" + differenceLoc + "\t" + modMatch.getString("peptideSequence").charAt(differenceLoc) + "->" + snpMatch.getString("peptideSequence").charAt(differenceLoc));
						break;	
					}
				}
			}
		}
	}
	
	/**
	 * We're checking that the two peptides are basically the same except they have one amino acid different. (from a SNP)
	 * @param matchA
	 * @param matchB
	 * @return returns -1 if they have more or less than one difference (eg no differences or more than one difference).  else returns the location of the difference
	 */
	private static int oneDifference(Match matchA, Match matchB) {
		String matchAString = matchA.getString("peptideSequence");
		String matchBString = matchB.getString("peptideSequence");
		if (matchAString.length() != matchBString.length()) return -1;
		int differenceCount = 0;
		int differenceLoc = -1;
		for (int i = 0; i < matchAString.length(); i++) {
			if (matchAString.charAt(i) != matchBString.charAt(i)) {
				differenceCount++;
				if (differenceCount > 1) return -1;
				differenceLoc = i;
			}	
		}
		return differenceLoc;
	}

}
