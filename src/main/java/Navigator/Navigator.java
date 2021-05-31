package Navigator;

import Peppy.U;

import java.io.*;
import java.util.ArrayList;

public class Navigator {
	
	public static void main(String [] args) {
		String reportFoldersName = "reportFolders.txt";
		U.p("loading " + reportFoldersName);
		File reportFoldersFile = new File(reportFoldersName);
		
		File unifiedReportsDir = new File("unifiedReports");
		unifiedReportsDir.mkdirs();
		File acidVariationDir = new File(unifiedReportsDir, "acidVariations");
		acidVariationDir.mkdirs();
		File modReportDir = new File(unifiedReportsDir, "modificatinReport");
		modReportDir.mkdirs();
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(reportFoldersFile));
			
			ArrayList<File> reportFolders = new ArrayList<File>();
			
			String line = br.readLine();
			while (line != null) {
				File file = new File(line);
				if (file.exists()) {
					reportFolders.add(file);
				} else {
					U.p("Could not find file: " + line);
				}
				line = br.readLine();
			}
			
			/* a list of our BestMatches */
			ArrayList<BestMatches> bestMatchesArray = new ArrayList<BestMatches>();
			
			/* a list of folders to ignore from our results */
			ArrayList<String> direcotryTitlesToIgnore = new ArrayList<String>();
			direcotryTitlesToIgnore.add("varimod");
//			direcotryTitlesToIgnore.add("HG19");
			direcotryTitlesToIgnore.add("mouse");
			direcotryTitlesToIgnore.add("xeno");
			direcotryTitlesToIgnore.add("gencode");
			direcotryTitlesToIgnore.add("personal");
			direcotryTitlesToIgnore.add("germline");
			direcotryTitlesToIgnore.add("subject");
			
			for (File folder: reportFolders) {
				BestMatches bestMatches = new BestMatches(folder, ResultsCategory.DNA, direcotryTitlesToIgnore);
				bestMatchesArray.add(bestMatches);
				
				File [] subFolders = folder.listFiles();
				for (File subFolder: subFolders) {
					if (subFolder.getName().toLowerCase().indexOf("varimod") != -1) {
						File varimodFile = new File(subFolder, "report.txt");
						if (varimodFile.exists()) {
							/* AA substitution report */
							File aaReportFile = new File(acidVariationDir, folder.getName() + ".html");
							AASubstitutionReport aasr = new AASubstitutionReport(varimodFile, aaReportFile);
							
							/* modification report */
							File sampleModDir  = new File(modReportDir, folder.getName());
							sampleModDir.mkdir();
							ModificationReport.findStickyModifications(aasr.getMatches(), sampleModDir);
						}
					}
				}
				
				
			}
			
			BestMatches.createUnifiedSamplesReport(bestMatchesArray, "peptideSequence", unifiedReportsDir);
			
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}		
		
		U.p("done");

	}

}
