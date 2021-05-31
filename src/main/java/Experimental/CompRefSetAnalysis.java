package Experimental;

import Navigator.MatchRow;
import Peppy.U;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

/*
 * Goals:
 * 1) find how many unique sequences per level
 * 2) print unique sequences from hg19 onward
 * 3) find increase identified in hg19 from other levels
 * 4) find increase in supporting sequences
 * 5) find set increases for hg19 discoveries
 * 6) set increase where only secondary peptides found
 */
public class CompRefSetAnalysis {
	static String uniprot;
	static Hashtable<String, String> inUniprot = new Hashtable<String, String>();
	static Hashtable<String, String> notInUniprot = new Hashtable<String, String>();
	static boolean filter = true;
	
	public static void main(String args[]) {
		init();
		createReport();
		U.p("done");
	}
	
	public static void init() {
		uniprot = U.readFileToString("/Users/risk2/PeppyData/public/sequences/protein/UniProt-HUMAN-20130918.fasta");
		uniprot = uniprot.replace("\r", "");
	}
	
	public static void createReport() {
		
		//load hg19 sequences from first analysis
		ArrayList<MatchRow> firstMatches = loadMatches("/Volumes/Research/CPTAC-CompRef/reports", "HG19");
		U.p("firstMatches " + firstMatches.size());
		Hashtable<String, MatchRow> firstHash = getMatchHash(firstMatches, "peptideSequence", null);
		U.p("firstHash " + firstHash.size());

//		MatchRow LETEIEALKEELLFMK = firstHash.get("LETEIEALKEELLFMK");
//		U.p(LETEIEALKEELLFMK.getInt("start") + "\t" + LETEIEALKEELLFMK.getInt("stop") + "\t" + LETEIEALKEELLFMK.getString("sequenceName"));
		
		//load secondary sequences
		ArrayList<MatchRow> secondMatches = loadMatches("/Volumes/Research/CPTAC-CompRef/reports-second", "sixFrameSequences");
		U.p("secondMatches " + secondMatches.size());
		Hashtable<String, MatchRow> secondMatchesHash = new Hashtable<String, MatchRow>();
		ArrayList<MatchRow> secondMatchesTrimmed = new ArrayList<MatchRow>();
		for (MatchRow match: secondMatches) {
			String peptideSequence = match.getString("peptideSequence");
			if (firstHash.get(peptideSequence) == null) {
				secondMatchesHash.put(peptideSequence, match);
				secondMatchesTrimmed.add(match);
			}
		}
		U.p("secondMatchesHash " + secondMatchesHash.size());
		
		Hashtable<String, MatchRow> secondMatchesTrimmedHash = getMatchHash(secondMatchesTrimmed, "peptideSequence", null);
		U.p("secondMatchesTrimmed " + secondMatchesTrimmed.size());
		U.p("secondMatchesTrimmedHash " + secondMatchesTrimmedHash.size());
		
		Hashtable<String, Hashtable<String, String>> firstSetHash = getSetHash(firstMatches);
		Hashtable<String, Hashtable<String, String>> secondSetHash = getSetHash(secondMatches);
		
		//finding support sequences
		Hashtable<String, Hashtable<String, String>> supportHash = new Hashtable<String, Hashtable<String, String>>();
		for (MatchRow match: secondMatchesTrimmed) {
			String peptideSequence = match.getString("peptideSequence");
			
			//getting primary peptide
			String sequenceName = match.getString("sequenceName");
			String [] sequenceChunks = sequenceName.split(";");
			String primary = sequenceChunks[sequenceChunks.length - 1].trim();
			primary = primary.substring("acidSequence:".length());
			
			//retrieving hash of supports
			Hashtable<String, String> supports = supportHash.get(primary);
			if (supports == null) supports = new Hashtable<String, String>();
			
			//adding ours in
			supports.put(peptideSequence, peptideSequence);
			supportHash.put(primary, supports);
			
			
		}
		
		//finding the proportion more spectra were identified from first stage to second stage
		Hashtable<String, Integer> firstStageTallies = new Hashtable<String, Integer>();
		Hashtable<String, Integer> secondStageTallies = new Hashtable<String, Integer>();
		for (MatchRow match: firstMatches) {
			String peptideSequence = match.getString("peptideSequence");
			Integer tally = firstStageTallies.get(peptideSequence);
			if (tally == null) {
				firstStageTallies.put(peptideSequence, 1);
			} else {
				firstStageTallies.put(peptideSequence, tally + 1);
			}
		}
		for (MatchRow match: secondMatches) {
			String peptideSequence = match.getString("peptideSequence");
			Integer tally = secondStageTallies.get(peptideSequence);
			if (tally == null) {
				secondStageTallies.put(peptideSequence, 1);
			} else {
				secondStageTallies.put(peptideSequence, tally + 1);
			}
		}
		Hashtable<Double, Integer> spectrumIncreases = new Hashtable<Double, Integer>();
		for (String peptide: firstStageTallies.keySet()) {
			int first = firstStageTallies.get(peptide);
			int second = secondStageTallies.get(peptide);
			double ratio = Math.ceil((double) second / first);
			Integer tally = spectrumIncreases.get(ratio);
			if (tally == null) {
				spectrumIncreases.put(ratio, 1);
			} else {
				spectrumIncreases.put(ratio, tally + 1);
			}
		}
		ArrayList<Double> ratios = new ArrayList<Double>(spectrumIncreases.keySet());
		Collections.sort(ratios);
		for (double ratio: ratios) {
			U.p(ratio + "\t" + spectrumIncreases.get(ratio));
		}
		PrintWriter spectrumIncreaseWriter;
		try {
			spectrumIncreaseWriter = new PrintWriter(new FileWriter("spectrum increases.txt"));
			for (String peptide: firstStageTallies.keySet()) {
				int first = firstStageTallies.get(peptide);
				int second = secondStageTallies.get(peptide);
				double ratio = (double) second / first;
				spectrumIncreaseWriter.println(peptide + "\t" + first + "\t" + second + "\t" + ratio);
			}
			spectrumIncreaseWriter.flush();
			spectrumIncreaseWriter.close();
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		
		
		
		
		//secondary set counts based on the peptide in the targeted sequence name
		Hashtable<String, Hashtable<String, String>> secondSetAll = new Hashtable<String, Hashtable<String, String>>();
		Hashtable<String, String> unambiguousSecond = new Hashtable<String, String>();
		for (MatchRow match: secondMatches) {
			String peptideSequence = match.getString("peptideSequence");
			String path = match.getFile("FilePath").getAbsolutePath();
			String [] pathChunks = path.split("/");
			String set = pathChunks[pathChunks.length - 2];
			
			//getting primary peptide
			String sequenceName = match.getString("sequenceName");
			String [] sequenceChunks = sequenceName.split(";");
			String primary = sequenceChunks[sequenceChunks.length - 1].trim();
			primary = primary.substring("acidSequence:".length());
//			if (firstHash.get(primary) == null) {
//				U.p(peptideSequence + "\t" + sequenceName);
//				unambiguousSecond.put(peptideSequence, peptideSequence);
//			}
			if (firstHash.get(peptideSequence) != null) primary = peptideSequence;
			
			Hashtable<String, String> setHash = secondSetAll.get(primary);
			if (setHash == null) setHash = new Hashtable<String, String> ();
			setHash.put(set, set);
			secondSetAll.put(primary, setHash);
		}
		
//		for (String pep: unambiguousSecond.keySet()) {
//			U.p(">" + pep);
//			U.p(pep);
//		}
		
		try {
			PrintWriter pw = new PrintWriter(new FileWriter("set increases.txt"));
			
			Hashtable<Double, Integer> increaseCounts = new Hashtable<Double, Integer>();
			for(String key: firstSetHash.keySet()) {
				int firstSize = firstSetHash.get(key).size();
				
				int secondSize = 0;
				if (secondSetHash.get(key) != null) secondSize = secondSetHash.get(key).size();
				
				MatchRow match = firstHash.get(key);
				
				//number of supports found
				int supportSize = 0;
				Hashtable<String, String> supports = supportHash.get(key);
				if (supports!= null) supportSize = supports.size();
				
				//number of sets, maybe some just secondary
				int secondSetAllCount = 0;
				Hashtable<String, String> secondSetAllHash = secondSetAll.get(key);
				if (secondSetAllHash != null) secondSetAllCount = secondSetAllHash.size();
				
				double ratio = (double) secondSize / firstSize;
				
				double ceiling = Math.ceil(ratio);
				Integer increaseCount = increaseCounts.get(ceiling);
				if (increaseCount == null) {
					increaseCounts.put(ceiling, 1);
				} else {
					increaseCounts.put(ceiling, increaseCount + 1);
				}
				
				if (secondSetAllCount == 0) {
					U.p("key: " + key);
//					for (MatchRow m: secondMatches) {
//						if (m.getString("peptideSequence").equals("ATPGHTGCLSPGCPDQPAR")) {
//							U.p(m.getFile("FilePath").getAbsolutePath());
//							U.p(m.getString("sequenceName"));
//						}
//						
//					}
//					return;
				}
				
				
				pw.println(key + "\t" + firstSize + "\t" + secondSize + "\t" + secondSetAllCount  + "\t" + ratio + "\t" + match.getString("sequenceName") + "\t" + match.getInt("start") + "\t" + supportSize);
			}
			
			pw.flush();
			pw.close();
			
			//printing the histogram of increase ratios
//			ArrayList<Double> ratios = new ArrayList<Double>(increaseCounts.keySet());
//			Collections.sort(ratios);
//			for (double ratio: ratios) {
//				U.p(ratio + "\t" + increaseCounts.get(ratio));
//			}
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}
	
	public static Hashtable<String, Hashtable<String, String>> getSetHash(ArrayList<MatchRow> matches) {
		Hashtable<String, Hashtable<String, String>> hash = new Hashtable<String, Hashtable<String, String>>();
		for (MatchRow match: matches) {
			String peptideSequence = match.getString("peptideSequence");
			String path = match.getFile("FilePath").getAbsolutePath();
			String [] pathChunks = path.split("/");
			String set = pathChunks[pathChunks.length - 2];
			Hashtable<String, String> setHash = hash.get(peptideSequence);
			if (setHash == null) {
				setHash = new Hashtable<String, String> ();
				hash.put(peptideSequence, setHash);
			}
			setHash.put(set, set);
		}
		return hash;
	}
	
	public static Hashtable<String, MatchRow> getMatchHash(ArrayList<MatchRow> matches, String matchElement, Hashtable<String, MatchRow> sequencesToIgnore) {
		Hashtable<String, MatchRow> hash = new Hashtable<String, MatchRow>();
		for (MatchRow match: matches) {
			String peptideSequence = match.getString(matchElement);
			if (sequencesToIgnore == null) {
				hash.put(peptideSequence, match);
			} else {
				if (sequencesToIgnore.get(peptideSequence) == null) {
					hash.put(peptideSequence, match);
				}
			}
		}
		return hash;
	}
	
	public static ArrayList<MatchRow> loadMatches(String filePath, String folderPattern) {
		File compref = new File(filePath);
		File [] reports = compref.listFiles();
		ArrayList<MatchRow> out = new ArrayList<MatchRow>();
		for (File reportFolder: reports) {
			if (reportFolder.isFile()) continue;
			if (reportFolder.isHidden()) continue;
			
			//find sub-folder that matches pattern, load that report
			File [] subFolders = reportFolder.listFiles();
			for (File subFolder: subFolders) {
				if (subFolder.getName().indexOf(folderPattern) != -1) {
					File reportFile = new File(subFolder, "report.txt");
					ArrayList<MatchRow> matches =  MatchRow.loadMatches(reportFile);
					if (filter) {
						for (MatchRow match: matches) {
							String peptideSequence = match.getString("peptideSequence");
							if (notInUniprot.get(peptideSequence) != null) {
								out.add(match);
							} else {
								if (inUniprot.get(peptideSequence) == null) {
									if (uniprot.indexOf(peptideSequence) != -1) {
										inUniprot.put(peptideSequence, peptideSequence);
									} else {
										out.add(match);
										notInUniprot.put(peptideSequence, peptideSequence);
									}
								}
							}
						}
					} else {
						out.addAll(matches);
					}
					break;
				}
			}

			
		}
		return out;
	}
	
	public static void getPostGenomeSequences() {
		File compref = new File("/Volumes/Research/CPTAC-CompRef/reports-second");
		File [] reports = compref.listFiles();
		
		Hashtable<String, Integer> supportSequenceCounts = new Hashtable<String, Integer>();
		Hashtable<String, String> supportSequenceRegions = new Hashtable<String, String>();
		
		for (File reportFolder: reports) {
			if (reportFolder.isFile()) continue;
			if (reportFolder.isHidden()) continue;
			
			File reportFile = new File(reportFolder, "1 sixFrameSequences/report.txt");
			ArrayList<MatchRow> matches =  MatchRow.loadMatches(reportFile);
			for (MatchRow match: matches) {
				String sequence = match.getString("peptideSequence");
				String sequenceName = match.getString("SequenceName");
				
				//Getting components of the name that tell about genomic location
				String [] nameChunks = sequenceName.split(";");
				String chr = nameChunks[0].trim();
				String start = nameChunks[3].split(":")[1].trim();
				String region = chr + start.substring(0,start.length() - 5);
				
				supportSequenceRegions.put(region, sequence);

				if (supportSequenceCounts.get(sequence) == null) {
					supportSequenceCounts.put(sequence,1);
				} else {
					int count = supportSequenceCounts.get(sequence);
					count++;
					supportSequenceCounts.put(sequence,count);
				}

				
			}
		}
		
		//print the results
		for (String region: supportSequenceRegions.keySet()) {
			String sequence = supportSequenceRegions.get(region);
			U.p(region + "\t" + sequence);
		}
	}
	
	
	public static Hashtable<String, Integer> getGenomeSequenceCounts() {

		
		File compref = new File("/Volumes/Research/CPTAC-CompRef/reports");
		File [] reports = compref.listFiles();
		
		Hashtable<String, Integer> proteinSequenceCounts = new Hashtable<String, Integer>();
		Hashtable<String, Integer> genomeSequenceCounts = new Hashtable<String, Integer>();
		
		
		for (File reportFolder: reports) {
			if (reportFolder.isFile()) continue;
			if (reportFolder.isHidden()) continue;
			
			File [] folders = reportFolder.listFiles();
			ArrayList<File> theSubFolders = new ArrayList<File>();
			for (int index = 0; index < folders.length; index++) {
				File folder = folders[index];
				if (folder.getName().startsWith("10")) continue;
				if (folder.getName().startsWith("11")) continue;
				if (folder.getName().startsWith("12")) continue;
				theSubFolders.add(folder);
			}
			for (int index = 0; index < folders.length; index++) {
				File folder = folders[index];
				if (folder.getName().startsWith("10")) theSubFolders.add(folder);
				if (folder.getName().startsWith("11")) theSubFolders.add(folder);
				if (folder.getName().startsWith("12")) theSubFolders.add(folder);
			}
			
			boolean beforeHG19 = true;
			boolean inHG19 = false;
			boolean pastHG19 = false;
			for (File subFolder: theSubFolders) {
				
				//skip if not a folder
				if (subFolder.isFile()) continue;
				
				//skip if no report
				File reportFile = new File(subFolder, "report.txt");
				if (!reportFile.exists()) continue;
				
				//if (subFolder.getName().indexOf("p5") != -1) continue;
				//if (subFolder.getName().indexOf("p6") != -1) continue;
				if (subFolder.getName().indexOf("varimod") != -1) continue;
				
				if (subFolder.getName().indexOf("HG19") != -1) inHG19 = true;
				if (inHG19 || pastHG19) beforeHG19 = false;
				if (beforeHG19) continue;
				
				ArrayList<MatchRow> matches =  MatchRow.loadMatches(reportFile);
				for (MatchRow match: matches) {
					String sequence = match.getString("peptideSequence");
					
					if (beforeHG19) {
						if (proteinSequenceCounts.get(sequence) == null) {
							proteinSequenceCounts.put(sequence,1);
						} else {
							int count = proteinSequenceCounts.get(sequence);
							count++;
							proteinSequenceCounts.put(sequence,count);
						}
					}
					
					if (inHG19) {
						/*only count sequence if not found previously */
						if (proteinSequenceCounts.get(sequence) == null) {
							if (genomeSequenceCounts.get(sequence) == null) {
								genomeSequenceCounts.put(sequence,1);
							} else {
								int count = genomeSequenceCounts.get(sequence);
								count++;
								genomeSequenceCounts.put(sequence,count);
							}
						}
					}
				}
				
				//if this was hg19, then we are now past it
				if (inHG19) {
					pastHG19 = true;
					inHG19 = false;
				}
				
				
			}
		}
		return genomeSequenceCounts;

	}
	
	/*
	 * 1) find how many unique sequences per level
	 * 2) print unique sequences from hg19 onward
	 */
	public static void UniquePerLevel() {
		
		
		try {
			PrintWriter sequenceCountWriter = new PrintWriter(new FileWriter("CPTAC sequence count report.txt"));
			PrintWriter blastReportWriter = new PrintWriter(new FileWriter("CPTAC blast report.txt"));
			PrintWriter postGenomeWriter = new PrintWriter(new FileWriter("CPTAC post genome.txt"));
			
			
			File compref = new File("/Volumes/Research/CPTAC-CompRef/reports");
			File [] reports = compref.listFiles();
			
			//where we save the sequences found after hg19
			Hashtable<String, String> postGenomeSequences = new Hashtable<String, String>();
			
			for (File reportFolder: reports) {
				if (reportFolder.isFile()) continue;
				if (reportFolder.isHidden()) continue;
				sequenceCountWriter.println(reportFolder.getName());
				U.p(reportFolder.getName());
				Hashtable<String, String> sequences = new Hashtable<String, String>();
				
				File [] folders = reportFolder.listFiles();
				ArrayList<File> theSubFolders = new ArrayList<File>();
				for (int index = 0; index < folders.length; index++) {
					File folder = folders[index];
					if (folder.getName().startsWith("10")) continue;
					if (folder.getName().startsWith("11")) continue;
					if (folder.getName().startsWith("12")) continue;
					theSubFolders.add(folder);
				}
				for (int index = 0; index < folders.length; index++) {
					File folder = folders[index];
					if (folder.getName().startsWith("10")) theSubFolders.add(folder);
					if (folder.getName().startsWith("11")) theSubFolders.add(folder);
					if (folder.getName().startsWith("12")) theSubFolders.add(folder);
				}
				
				boolean pastHG19 = false;
				for (File subFolder: theSubFolders) {
					
					//skip if not a folder
					if (subFolder.isFile()) continue;
					
					//skip if no report
					File reportFile = new File(subFolder, "report.txt");
					if (!reportFile.exists()) continue;
					
					if (subFolder.getName().indexOf("p5") != -1) continue;
					if (subFolder.getName().indexOf("p6") != -1) continue;
					if (subFolder.getName().indexOf("varimod") != -1) continue;
					
					
					//count unique sequences
					int uniqueSequenceCount = 0;
					
					ArrayList<MatchRow> matches =  MatchRow.loadMatches(reportFile);
					for (MatchRow match: matches) {
						String sequence = match.getString("peptideSequence");
						if (sequences.get(sequence) == null) {
							sequences.put(sequence, sequence);
							uniqueSequenceCount++;
							
							if (pastHG19) {
								if (postGenomeSequences.get(sequence) == null) {
									postGenomeSequences.put(sequence, match.getString("sequenceName"));
								}
							}
						}
					}
					
					sequenceCountWriter.println(subFolder.getName() + "\t" + uniqueSequenceCount);
					
					//if this was hg19, then we are now past it
					if (subFolder.getName().indexOf("HG19") != -1) pastHG19 = true;
					
					
				}
				
				sequenceCountWriter.println();
			}
			
			
			//printing the sequences found after hg19
			for(String key: postGenomeSequences.keySet()) {
				blastReportWriter.println(">" + key);
				blastReportWriter.println(key);
				blastReportWriter.println();
				
				postGenomeWriter.println(key + "\t" + postGenomeSequences.get(key));
			}
			
			
			
			
			sequenceCountWriter.flush();
			sequenceCountWriter.close();
			blastReportWriter.flush();
			blastReportWriter.close();
			postGenomeWriter.flush();
			postGenomeWriter.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		U.p("done");
	}

}
