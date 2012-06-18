package Navigator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

import Math.MassError;
import Peppy.Definitions;
import Peppy.ModificationEntry;
import Peppy.Properties;
import Peppy.U;

/**
 * Copyright 2012, Brian Risk
 * Released under the Netscape Public License
 * 
 * @author Brian Risk
 *
 */
public class ModificationReport {
	
	
	public static void main(String args[]) {
		new ModificationReport();
		U.p("done");
	}
	
	public ModificationReport() {
		double matchScoreCutoff = 28.24;
		Properties.precursorTolerance = 25;
		
		/* load the sample */
		Sample whim16_41 = new Sample("WHIM16 - 41", Sample.REFERENCE_PROTEIN);
		whim16_41.loadResults(new File("//Users/risk2/PeppyData/WashU/reports/41/WHIM16/WashU WHIM16 human/report.txt"));
		whim16_41.reduceMatches(matchScoreCutoff);
		ArrayList<Modification> whim2Mods = getSampleModifications(whim16_41);
		
//		for (Modification mod: whim2Mods) {
//			U.p(mod.getScore());
//		}
		
		Sample whim16_33 = new Sample("WHIM16 - 33", Sample.REFERENCE_PROTEIN);
		whim16_33.loadResults(new File("//Users/risk2/PeppyData/WashU/reports/33/WHIM16/WashU WHIM16 human/report.txt"));
		whim16_33.reduceMatches(matchScoreCutoff);
		ArrayList<Modification> whim16Mods = getSampleModifications(whim16_33);
		
		
		/* the modification lists should have all the same elements and be in the same order
		 * so we're going to take advantage of that
		 */
		
		ArrayList<ModificationComparison> comparisons = new ArrayList<ModificationComparison>(whim2Mods.size());
		for (int i = 0; i < whim2Mods.size(); i++) {
			Modification mod2 = whim2Mods.get(i);
			Modification mod16 = whim16Mods.get(i);
			ModificationComparison comparison = new ModificationComparison("" + mod2.getMass()); 
			comparison.addModification(mod2);
			comparison.addModification(mod16);
			comparisons.add(comparison);
		}
//		Collections.sort(comparisons);
		
		/* find the general ratio.  Perhaps one analysis is more PTM friendly than another? */
		double scoreA = (double) whim2Mods.size() / whim16_41.getMatches().size();
		double scoreB = (double) whim16Mods.size() / whim16_33.getMatches().size();
		
		
		
		/* save our reports */
		File reportDir = new File("modification comparison");
		reportDir.mkdir();
		
		/* create main index file */
		try {

			PrintWriter matchWriter = new PrintWriter(new BufferedWriter(new FileWriter(new File(reportDir, "index.html"))));
			matchWriter.println(HTML.sortableTableHeader);
			matchWriter.println("The ratio of set A to B is: " + (scoreA / scoreB));
			matchWriter.println(HTML.tableTop);
			matchWriter.println("<tr>");
			matchWriter.println("<th>mass</th>");
			matchWriter.println("<th>score</th>");
			matchWriter.println("<th>dominant</th>");
			matchWriter.println("<th>WHIM16 - 41</th>");
			matchWriter.println("<th>WHIM16 - 33</th>");
			matchWriter.println("</tr>");
			matchWriter.println("</thead>");
			matchWriter.println("<tbody>");
			for (int i = 0; i < comparisons.size(); i++) {
				ModificationComparison comparison = comparisons.get(i);
				if (comparison.getMinMatchCount() < 10) continue;
//				if (comparison.getScore() < 1.3) continue;
				matchWriter.println("<tr>");
				matchWriter.println("<td>"  + comparison.getName() + "</td>");
				matchWriter.println("<td>" + comparison.getScore() + "</td>");
				matchWriter.println("<td>" + comparison.getDominantSampleName() + "</td>");
				
				Modification modification =  whim2Mods.get(i);
				String fileName = modification.getSample().getName() + modification.getMass() + ".html";
				File saveFile = new File(reportDir, fileName);
//				createSampleReport(saveFile, whim2Mods);
				matchWriter.println("<td>" + modification.getMatchesSize() + "</td>");
				
				modification =  whim16Mods.get(i);
				fileName = modification.getSample().getName() + modification.getMass() + ".html";
				saveFile = new File(reportDir, fileName);
//				createSampleReport(saveFile, whim16Mods);
				matchWriter.println("<td>" + modification.getMatchesSize() + "</td>");
				

	
				matchWriter.println("</tr>");
			}
			matchWriter.println(HTML.sortableTableFooter);
			matchWriter.flush();
				matchWriter.close();
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
		
		
	}
	
	private void createSampleReport(File file, ArrayList<Modification> modifications) {
		/* create main index file */
		try {

			PrintWriter matchWriter = new PrintWriter(new BufferedWriter(new FileWriter(file)));
			matchWriter.println(HTML.sortableTableHeader);
			matchWriter.println(HTML.tableTop);
			matchWriter.println("<tr>");
			matchWriter.println("<th>count</th>");
			matchWriter.println("<th>mass</th>");
			matchWriter.println("<th>possibilities</th>");
			matchWriter.println("</tr>");
			matchWriter.println("</thead>");
			matchWriter.println("<tbody>");
			
			for (Modification modification: modifications) {
				if (modification.getMatchesSize() > 0) {
					matchWriter.println("<tr>");
					matchWriter.println("<td>" + modification.getMatchesSize() + "</td>");
					matchWriter.println("<td>" + modification.getMass() + "</td>");
					
					/* printing out the modification descriptions */
					ArrayList<ModificationEntry> entries = modification.getEntries();
					matchWriter.println("<td><ul>");
					for (ModificationEntry entry: entries) {
						matchWriter.println("<li>" + entry.getDescription() + "</li>");
					}
					
					matchWriter.println("</ul></td>");
					matchWriter.println("</tr>");
				}
			}
			
			matchWriter.println(HTML.sortableTableFooter);
			matchWriter.flush();
			matchWriter.close();
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}
	
	
	private ArrayList<Modification> getSampleModifications(Sample sample) {
		/* create a set of known modifications */
		Hashtable<Double, Modification> modHash = new Hashtable<Double, Modification>();
		for (ModificationEntry mod: Definitions.modificationEntries) {
			Modification modification = modHash.get(mod.getMonoMass());
			if (modification == null) {
				modification = new Modification(mod.getMonoMass(), sample);
				modHash.put(mod.getMonoMass(), modification);
			}
			modification.addEntry(mod);
		}
		ArrayList<Modification> modifications = new ArrayList<Modification>(modHash.values());
		Collections.sort(modifications);
		
		
		for (Match match: sample.getMatches()) {
			if (match.getBoolean("isModified")) {
				double peptideMass = match.getDouble("PrecursorNeutralMass");
				double modMass = match.getDouble("modMass");
				double error = MassError.getDaltonError(Properties.fragmentTolerance, peptideMass);
				int start = Math.abs(Collections.binarySearch(modifications, new Modification(modMass - error, sample))) - 1;
				int stop = Math.abs(Collections.binarySearch(modifications, new Modification(modMass + error, sample))) - 1;
				for (int i = start; i < stop; i++) {
					Modification modification = modifications.get(i);
					modification.addMatch(match);
					modification.addToScore(1000000.0 / sample.getMatches().size());
				}

			}
		}
		return modifications;
	}
	
	

}
