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
import Reports.HTMLPage;

/**
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class ModificationReport {
	
	
	public static void main(String args[]) {
		
		
		
		ArrayList<Match> matches = Match.loadMatches(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/RWA1/4 RWA1 - varimod/report.txt"));
		
		
		
//		ArrayList<Match> matches = new ArrayList<Match>();
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/Documents/MD Anderson/Visit results/group1/3 group1 - varimod/report.txt")));
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/Documents/MD Anderson/Visit results/group2/3 group2 - varimod/report.txt")));
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/Documents/MD Anderson/Visit results/group3/3 group3 - varimod/report.txt")));
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/Documents/MD Anderson/Visit results/group4/3 group4 - varimod/report.txt")));
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/Documents/MD Anderson/Visit results/group5/3 group5 - varimod/report.txt")));
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/Documents/MD Anderson/Visit results/group6/3 group6 - varimod/report.txt")));

		
		
//		ArrayList<Match> matches = new ArrayList<Match>();
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM2-Ellis033/5 WHIM2 - varimod/report.txt")));
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM2-Ellis041/6 WHIM2 - varimod/report.txt")));
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/PeppyData/CPTAC/reports/WHIM2-Ellis043/5 WHIM2 - varimod/report.txt")));
//		matches.addAll(Match.loadMatches(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM2-Ellis043/5 WHIM2 - varimod/report.txt")));
		
		File saveFolder = new File("modification report");
		findStickyModifications(matches , saveFolder);
		createModSpreadsheet(matches , saveFolder);
		createModificationTable(matches , saveFolder);
		
		U.p("done");
	}
	
	public static void createModSpreadsheet(File reportFile, File saveFolder) {
		ArrayList<Match> matches = Match.loadMatches(reportFile);
		createModSpreadsheet(matches, saveFolder);
	}
	
	public static void createModSpreadsheet(ArrayList<Match> matches, File saveFolder) {

		
		try {
			
			for (int targetModMass = 0; targetModMass < 101; targetModMass++) {
				Hashtable<Character, Integer> acidModificationTallies = new Hashtable<Character, Integer>();
				int nTerminalCount = 0;
				int cTerminalCount = 0;
				
				
				
				for (Match match: matches) {
					if (match.getBoolean("isModified")  && match.getBoolean("modLocCertain") ) {
						int modMass = (int) Math.round(match.getDouble("modMass"));
						if (modMass == targetModMass) {
							String peptideSequence = match.getString("peptideSequence");
							int modIndex = match.getInt("modIndex");
							
							/*
							 * skip if n-terminal
							 */
							//if (modIndex == 0) continue;
							
							char acid = peptideSequence.charAt(modIndex);
							
							Integer acidTally = acidModificationTallies.get(acid);
							if (acidTally == null) {
								acidModificationTallies.put(acid, 1);
							} else {
								acidModificationTallies.put(acid, acidTally + 1);
							}
							
							if (modIndex == 0) nTerminalCount++;
							if (modIndex == peptideSequence.length() - 1) cTerminalCount++;
						}
						
					}
				}
				
				
				/* find if there is a prevailing residue location */
				ArrayList<Integer> values = new ArrayList<Integer>(acidModificationTallies.values());
				
				int sum = 0; 
				for (Integer value: values) {
					sum += value;
				}
				if (sum < 1) continue;
				int max = Collections.max(values);
				double maxRatio = (double) max / sum ;
				
//				if (maxRatio < .2) continue;
				
//				File saveFolder = new File(reportFile.getParentFile().getParentFile().getName() + " mods");
				saveFolder.mkdir();
//				PrintWriter pw = new PrintWriter(new FileWriter(new File(saveFolder, targetModMass  + " (" + sum + ") modAcids.txt")));
				PrintWriter pw = new PrintWriter(new FileWriter(new File(saveFolder, sum  + " (" + targetModMass + ") modAcids.txt")));
				
				/*
				 * Print the amino acids with their modification counts
				 */
				ArrayList<Character> keys = new ArrayList<Character>(acidModificationTallies.keySet());
				Collections.sort(keys);
				for (Character key: keys) {
					int tally =  acidModificationTallies.get(key);
					double ratio = (double) tally / sum;
					pw.println(key + "\t" + tally  + "\t" + ratio);
				}
				
				/*
				 * print the N and C terminal counts
				 */
				pw.println();
				pw.println("n-terminal count: " + nTerminalCount + "(" + ((double) nTerminalCount / sum) + ")");
				pw.println("c-terminal count: " + cTerminalCount + "(" + ((double) cTerminalCount / sum) + ")");
				
				/*
				 * print the modifications this might be
				 */
				pw.println();
				for (ModificationEntry mod: Definitions.modificationEntries) {
					if (Math.round(mod.getMonoMass()) == targetModMass) {
						pw.println(mod.getDescription());
					}
				}
				
				pw.flush();
				pw.close();
				
				
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	
	public static void createModificationTable(ArrayList<Match> matches, File saveFolder) {
		
		/*
		 * First we must get amino acid frequencies.  This will help us know the
		 * probability of a modification randomly landing on a given acid.
		 */
		int totalAminoAcids = 0;
		Hashtable<Character, Integer> acidCounts = new Hashtable<Character, Integer>();
		
		for (Match match: matches) {
			String peptideSequence = match.getString("peptideSequence");
			char [] peptideAcids = peptideSequence.toCharArray();
			
			totalAminoAcids += peptideSequence.length();
			
			/* add in all acids except the n-terminal*/
			char acid;
			for (int i = 1; i < peptideAcids.length; i++) {
				acid = peptideAcids[i];
				/* ignore stops */
				if (acid == '.') continue;
				
				Integer acidCount = acidCounts.get(acid);
				if (acidCount == null) acidCount = 0;
				acidCount++;
				acidCounts.put(acid, acidCount);
			}
		}
		
		/*
		 * now to find the modification frequencies
		 */
		Hashtable<Integer, Integer> massTallies = new Hashtable<Integer, Integer>();
		Hashtable<String, ModificationType> modificationTypes = new Hashtable<String, ModificationType>();
		for (Match match: matches) {
			if (match.getBoolean("isModified")  && match.getBoolean("modLocCertain") ) {
//			if (match.getBoolean("isModified") ) {
				int modIndex = match.getInt("modIndex");
				
//				if (modIndex == 0) continue; 
				int modMass = (int) Math.round(match.getDouble("modMass"));
				String peptideSequence = match.getString("peptideSequence");
				char acid = peptideSequence.charAt(modIndex);
				
				
				
				/* keeping track of how many times we have seen a modification of this mass */
				Integer massTally = massTallies.get(modMass);
				if (massTally == null) massTally = 0;
				massTally++;
				massTallies.put(modMass, massTally);
				
				/* tracking how many mods of this mass at the given acid were seen */
				String hashString = acid + "" + modMass;
				ModificationType modificationType = modificationTypes.get(hashString);
				if (modificationType == null) modificationType = new ModificationType(acid, modMass);
				modificationType.incrementTally();
				modificationTypes.put(hashString, modificationType);
				
			}
				
		}
		
		/*
		 * calculate the scores
		 */
		ArrayList<ModificationType> modificationTypeArray = new ArrayList<ModificationType>(modificationTypes.values());
		for (ModificationType modificaitonType: modificationTypeArray) {
			int totalModOccurrences = massTallies.get(modificaitonType.getMass());
			double acidPercentage = (double) acidCounts.get(modificaitonType.getAcid()) / totalAminoAcids;
			modificaitonType.calculateScore(totalModOccurrences, acidPercentage);
		}
		
		Collections.sort(modificationTypeArray);
		
		/*
		 * print results
		 */
		saveFolder.mkdirs();
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(saveFolder,"mod score report.txt"))));
			for (ModificationType modificaitonType: modificationTypeArray) {
				pw.println(modificaitonType.getAcid() + "\t" + modificaitonType.getMass() + "\t" + modificaitonType.getScore()+ "\t" + modificaitonType.getTally());
			}
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		

		
	}
	
	
	
	public static void findStickyModifications(ArrayList<Match> matches, File saveFolder) {
		
		/*
		 * Want to count how many peptides contain a given acid.  Counting
		 * only once as blind modification only accounts for one acid
		 * being modified.  Counting all the occurrences of an acid will throw off
		 * the probabilities we are calculating.
		 */
		int totalAminoAcids = 0;
		Hashtable<Character, Integer> acidCounts = new Hashtable<Character, Integer>();
		
		for (Match match: matches) {
			String peptideSequence = match.getString("peptideSequence");
			char [] peptideAcids = peptideSequence.toCharArray();
			
			totalAminoAcids += peptideSequence.length();
			Hashtable<Character, Character> peptideHashtable = new Hashtable<Character, Character>();
			for (char acid: peptideAcids) {
				/* ignore stops */
				if (acid == '.') continue;
				peptideHashtable.put(acid, acid);
			}
			
			for (char acid: peptideHashtable.values()) {
				Integer acidCount = acidCounts.get(acid);
				if (acidCount == null) acidCount = 0;
				acidCount++;
				acidCounts.put(acid, acidCount);
			}
		}
		
		
		/*
		 * now to find the modification frequencies
		 */
		Hashtable<Integer, Integer> massTallies = new Hashtable<Integer, Integer>();
		Hashtable<String, ModificationEvent> modificationEvents = new Hashtable<String, ModificationEvent>();
		Hashtable<String, ModificationType> modificationTypes = new Hashtable<String, ModificationType>();
		Hashtable<String, Integer> sequenceModificaitonCounts = new Hashtable<String, Integer>();
		for (Match match: matches) {
			if (match.getBoolean("isModified")  && match.getBoolean("modLocCertain") ) {
				int modIndex = match.getInt("modIndex");
//				if (modIndex == 0) continue; /* skip if n-terminal */
				int modMass = (int) Math.round(match.getDouble("modMass"));
				String peptideSequence = match.getString("peptideSequence");
				char acid = peptideSequence.charAt(modIndex);
				if (acid == '.') continue;
				
				
				
				/* keeping track of how many times we have seen a modification of this mass */
				Integer massTally = massTallies.get(modMass);
				if (massTally == null) massTally = 0;
				massTally++;
				massTallies.put(modMass, massTally);
				
				/* tracking how many mods of this mass at the given acid were seen */
				String hashString = peptideSequence + "-" + modIndex + "-" + modMass;
				ModificationEvent modificationEvent = modificationEvents.get(hashString);
				if (modificationEvent == null) {
					modificationEvent = new ModificationEvent(peptideSequence, modIndex, modMass, acid, match);
					
					/*if this is the first time we're seeing a mod event, incremnt that the sequence was modified */
					Integer sequenceModCount = sequenceModificaitonCounts.get(match.getString("SequenceName"));
					if (sequenceModCount == null) sequenceModCount = 0;
					sequenceModCount++;
					sequenceModificaitonCounts.put(match.getString("SequenceName"), sequenceModCount);
				}
				modificationEvent.incrementCount();
				modificationEvents.put(hashString, modificationEvent);
				
				/* tracking how many mods of this mass at the given acid were seen */
				String typeHash = acid + "" + modMass;
				ModificationType modificationType = modificationTypes.get(typeHash);
				if (modificationType == null) modificationType = new ModificationType(acid, modMass);
				modificationType.incrementTally();
				modificationTypes.put(typeHash, modificationType);
				
			}
				
		}
		
		/*
		 * find statistics:
		 * the mean for each modification type
		 */
		for (ModificationEvent modEvent: modificationEvents.values()) {
			String typeHash = modEvent.getAcid() + "" + modEvent.getModMass();
			ModificationType modificationType = modificationTypes.get(typeHash);
			int peptideWithAcidCount = acidCounts.get(modEvent.getAcid());
			double probability = (double) modificationType.getTally() / peptideWithAcidCount;
			modEvent.calculateScore(probability);
		}
		
		/*
		 * list of substitutions and their masses
		 */
		ArrayList<AASubstitution> substitutions = AASubstitution.generateListOfAASubstitutions();
//		ArrayList<AASubstitution> substitutions = AASubstitution.generateListOfSingleNucleotideAASubstitutions();
		
		Hashtable<String, AASubstitution> substitutionHashtable = new Hashtable<String, AASubstitution> (substitutions.size());
		for (AASubstitution substitution: substitutions) {
			String subHash = "" + substitution.getPrevious() + (int) Math.round(substitution.getDifference());
			substitutionHashtable.put(subHash, substitution);
		}
		
		
		
		ArrayList<ModificationEvent> modEventErray = new ArrayList<ModificationEvent>(modificationEvents.values()); 
		Collections.sort(modEventErray);
		
		/*
		 * print results
		 */
		saveFolder.mkdirs();
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(saveFolder,"sticky modification report.html"))));
			pw.println("<html><body><table border=1>");
			pw.println("<tr>");
			pw.println("<th>sequence</th>");
			pw.println("<th>seqModCount</th>");
			pw.println("<th>peptide</th>");
			pw.println("<th>acid</th>");
			pw.println("<th>new acid</th>");
			pw.println("<th>location</th>");
			pw.println("<th>mass shift</th>");
			pw.println("<th>count</th>");
			pw.println("<th>score</th>");
			for (ModificationEvent modificationEvent: modEventErray) {
				/* only want matches of a certain score threshold */
//				if (modificationEvent.getScore() < 20) break;
				
				/* trim the sequence name if it has a pipe character */
				String sequenceName = modificationEvent.getMatch().getString("SequenceName") ;
				if (sequenceName.indexOf("|") != -1) {
					sequenceName = sequenceName.split("\\|")[0];
				}
				
				/* seeing if it is an acid substitution */
				String subHash = "" + modificationEvent.getAcid() + modificationEvent.getModMass();
				AASubstitution substitution = substitutionHashtable.get(subHash);
				if (substitution == null) continue;
				
				pw.println("<tr>");
				
				pw.println("<td><a href=\"http://www.genecards.org/cgi-bin/carddisp.pl?gene=" + sequenceName +  "\">" + sequenceName + "</a></td>");
//				pw.println("<td><a href=\"http://www.uniprot.org/uniprot/" + sequenceName +  "\">" + sequenceName + "</a></td>");
				pw.println("<td>" + sequenceModificaitonCounts.get(modificationEvent.getMatch().getString("SequenceName")) + "</td>");
				pw.println("<td>" + modificationEvent.getPeptideSequence() + "</td>");
				pw.println("<td>" + modificationEvent.getAcid() + "</td>");
				pw.println("<td>" + substitution.getPresent() + "</td>");
				pw.println("<td>" + modificationEvent.getLocation() + "</td>");
				pw.println("<td>" + modificationEvent.getModMass() + "</td>");
				pw.println("<td>" + modificationEvent.getCount() + "</td>");
				pw.println("<td>" + modificationEvent.getScore() + "</td>");
				
			}
			pw.println("</table></body></html>");
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		

		
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
