package MoonShot;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;

import Math.MathFunctions;
import Navigator.Match;
import Navigator.Region;
import Peppy.Definitions;
import Peppy.Peptide;
import Peppy.Protein;
import Peppy.Sequence_Protein;
import Peppy.U;

public class MoonShot {
	
	public static void main(String [] args) {
		createRegionReport();
	}
	
	public static void createRegionReport() {
		
		File reportsFolder;
		String novelPath;
		boolean getGencode, getNovelDistance;
		
//		reportsFolder = new File("reports");
//		novelPath = "2 HG19";
//		getGencode = true;
//		getNovelDistance = true;
//		ArrayList<Match> firstPassMatches = getUnifiedMatches(reportsFolder, novelPath, getGencode, getNovelDistance);
//		printMatches(firstPassMatches, "master moon matches.txt");
//		
//		reportsFolder = new File("reports-complete/moon Shot second pass");
//		novelPath = "2 lung28";
//		getGencode = true;
//		getNovelDistance = true;
//		ArrayList<Match> targetedMatches = getUnifiedMatches(reportsFolder, novelPath, getGencode, getNovelDistance);
//		printMatches(targetedMatches, "targeted matches.txt");
		
//		reportsFolder = new File("MoonShot second pass");
//		novelPath = "2 01";
//		getGencode = true;
//		getNovelDistance = true;
//		ArrayList<Match> targetedMatches = getUnifiedMatches(reportsFolder, novelPath, getGencode, getNovelDistance);
//		printMatches(targetedMatches, "targeted matches.txt");
		
		String studyName = "CPTAC compRef";
		
		reportsFolder = new File("reports");
		novelPath = "HG19";
		getGencode = true;
		getNovelDistance = true;
		ArrayList<Match> targetedMatches = getUnifiedMatches(reportsFolder, novelPath, getGencode, getNovelDistance);
		printMatches(targetedMatches, studyName + " targeted matches.txt");
		
//		reportsFolder = new File("reports");
//		novelPath = "2 kentsisRegions";
//		getGencode = true;
//		getNovelDistance = true;
//		ArrayList<Match> targetedMatches = getUnifiedMatches(reportsFolder, novelPath, getGencode, getNovelDistance);
//		printMatches(targetedMatches, "kentsis targeted matches.txt");
		
		
//		ArrayList<Match> firstPassMatches = Match.loadMatches(new File("master moon matches.txt"));
//		assignGencodeAnnotations(firstPassMatches);
		
//		ArrayList<Match> targetedMatches =  Match.loadMatches(new File(studyName + " targeted matches.txt"));
		
		
		
//		remapTargetedMatches(targetedMatches);
//		assignP(targetedMatches);
		assignGencodeAnnotations(targetedMatches);
//		printMatches(targetedMatches, studyName + " remapped matches.txt");
		
		ArrayList<Match> unifiedMatches =  new ArrayList<Match>();
		unifiedMatches.addAll(targetedMatches);
//		unifiedMatches.addAll(firstPassMatches);
		ArrayList<MoonShotRegion> regions = new ArrayList<MoonShotRegion>();
		addMatchesToRegions(unifiedMatches, regions);
		
		//assign IDs to regions
		for (int id = 0; id < regions.size(); id++) {
			regions.get(id).setId(id);
		}
		
		//find how regions overlap
		//see how regions overlap with members
		for (int indexA = 0; indexA < regions.size() - 1; indexA++) {
			MoonShotRegion regionA = regions.get(indexA);
			for (int indexB = indexA + 1; indexB < regions.size(); indexB++) {
				MoonShotRegion regionB = regions.get(indexB);
				regionA.getUniqueCount(regionB);
				regionB.getUniqueCount(regionA);
			}
		}
		
		
		
		//save regions to a file
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(studyName + " regions.txt")));
			for (MoonShotRegion region: regions) {
				pw.println(region);
			}
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//save html
		File reportDir = new File(studyName + " regions report");
		File regionsDir = new File (reportDir, "regions");
		regionsDir.mkdirs();
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(reportDir, "index.html"))));
			pw.println("<html><body><table>");
			pw.println("<thead><tr>");
			pw.println("<th>ID</th>");
			pw.println("<th>Gene</th>");
			pw.println("<th>Note</th>");
			pw.println("<th>Location</th>");
			pw.println("<th>Matches</th>");
			pw.println("<th>Sequences</th>");
			pw.println("<th>Distance</th>");
			pw.println("<th>Samples</th>");
			pw.println("<th>Uniquet</th>");
			pw.println("<th>Similar</th>");
			pw.println("</tr></thead>");
			
			for (MoonShotRegion region: regions) {
				pw.println(region.toHTML());
			}
			pw.println("</body></html>");
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(reportDir, "bed.txt"))));
			
			for (Match match: unifiedMatches) {
				pw.println(match.toBED());
			}
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
		
		
		U.p("done");
	}
	
	public static void assignP(ArrayList<Match> matches) {
		Hashtable<String, String> nonSmokers = new Hashtable<String, String>();
		
		
		nonSmokers.put("H1299", "xxx");
		nonSmokers.put("H838", "xxx");
		nonSmokers.put("H1437", "xxx");
		nonSmokers.put("H1568", "xxx");
		nonSmokers.put("H1693", "xxx");
		nonSmokers.put("H1993", "xxx");
		nonSmokers.put("H522", "xxx");
		nonSmokers.put("H2009", "xxx");
		nonSmokers.put("H23", "xxx");
		nonSmokers.put("H1944", "xxx");
		nonSmokers.put("HCC4019", "xxx"); 
		nonSmokers.put("H1792", "xxx");
		nonSmokers.put("HCC4017", "xxx");
		nonSmokers.put("H1355", "xxx");



		
		ArrayList<String> cellLines = new ArrayList<String>();
		cellLines.add("H1651");
		cellLines.add("H1435");
		cellLines.add("H2342");
		cellLines.add("H1563");
		cellLines.add("H1838");
		cellLines.add("H1975");
		cellLines.add("H1650");
		cellLines.add("H3255");
		cellLines.add("HCC827");
		cellLines.add("HCC4011");
		cellLines.add("ipas339"); //HCC2279
		cellLines.add("H820");
		cellLines.add("H1299");
		cellLines.add("H838");
		cellLines.add("H1437");
		cellLines.add("H1568");
		cellLines.add("H1693");
		cellLines.add("H1993");
		cellLines.add("H522");
		cellLines.add("H2009");
		cellLines.add("H23");
		cellLines.add("H1944");
		cellLines.add("HCC4019");
		cellLines.add("H1792");
		cellLines.add("H1355");
		cellLines.add("HCC4017");
		
		final double P = (double) nonSmokers.size() / cellLines.size();
		
		//generate list of unique peptides
		Hashtable<String, ArrayList<Match>> peptides = new Hashtable<String, ArrayList<Match>>(); 
		for (Match match: matches) {
			String peptide = match.getString("peptideSequence");
			ArrayList<Match> peptideMatches = peptides.get(peptide);
			if (peptideMatches == null) peptideMatches = new ArrayList<Match>();
			peptideMatches.add(match);
			peptides.put(peptide, peptideMatches);
			if (peptide.equals("LGLLVNIYR")) U.p(peptide + ": " + match.getInt("start"));
			if (peptide.equals("GAYGGGYGGYDDYGGYNDGYGFGSDR")) U.p(peptide + ": " + match.getInt("start"));
			if (peptide.equals("DDCGDGSDEASCPK")) U.p(peptide + ": " + match.getInt("start"));
			if (peptide.equals("VLQQVLER")) U.p(peptide + ": " + match.getInt("start"));
			
		}
		
		/*
		 * for each peptide, we need the tally for:
		 * 1) the number of nonSmoker cell lines it is found in
		 * 2) the total number of cell lines it is found in
		 */
		Hashtable<String, Double> peptidePValues = new Hashtable<String, Double>();
		for (String peptide: peptides.keySet()){
			Hashtable<String, String> keyPeptides = new Hashtable<String, String>();
			Hashtable<String, String> cellLinePeptides = new Hashtable<String, String>();
			for (Match match: peptides.get(peptide)) {
				String fileName = match.getFile("FilePath").getAbsolutePath();
				for (String cellLine: cellLines) {
					if (fileName.toLowerCase().indexOf(cellLine.toLowerCase()) == -1) continue;
					cellLinePeptides.put(cellLine, cellLine);
					if (nonSmokers.get(cellLine) != null) {
						keyPeptides.put(cellLine, cellLine);
					}
				}
			}
			double peptideP = MathFunctions.approximateBinomialProbability(cellLinePeptides.size(), keyPeptides.size(), P);
			peptidePValues.put(peptide, peptideP);
			U.p(peptide + "\t" + peptideP + "\t" + cellLinePeptides.size() + "\t" + keyPeptides.size() );

				
		}
	}
	
	
	public static void tally() {
		ArrayList<String> nonSmokers = new ArrayList<String>();
		nonSmokers.add("H1651");
		nonSmokers.add("H1435");
		nonSmokers.add("H2342");
		nonSmokers.add("H1563");
		nonSmokers.add("H1838");
		nonSmokers.add("H1975");
		nonSmokers.add("H1650");
		nonSmokers.add("H3255");
		nonSmokers.add("HCC827");
		nonSmokers.add("HCC4011");
		nonSmokers.add("ipas339"); //HCC2279
		nonSmokers.add("H820");
		
		
		File reportsFolder = new File("MoonShot second pass");	
		File [] reportsFolders = reportsFolder.listFiles();
		int tally = 0;
		int nonSmokerTally = 0;
		for (File reportFolder: reportsFolders) {
			File metrics = new File(reportFolder, "metrics.txt");
			if (!metrics.exists()) continue;
			try {
				BufferedReader metricsReader = new BufferedReader(new FileReader(metrics));
				String line = metricsReader.readLine();
				int spectraIdentified = 0;
				while (line != null) {
					String [] chunks = line.split(" ");
					if (chunks[0].equals("spectra")) {
						spectraIdentified = Integer.parseInt(chunks[2]);
					}
					line = metricsReader.readLine();
				}
				tally += spectraIdentified;
				
				for (String non: nonSmokers) {
					if (reportFolder.getName().toLowerCase().indexOf(non.toLowerCase()) == -1) continue;
					nonSmokerTally += spectraIdentified;
				}
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		U.p(tally);
		U.p(nonSmokerTally);
	}
	
	public static void addMatchesToRegions(ArrayList<Match> matches, ArrayList<MoonShotRegion> regions) {
		for (Match match: matches) {
			
			//try adding to all of the existing regions
			boolean wasAdded = false;
			for (MoonShotRegion region: regions) {
				boolean added = region.addMatch(match);
				if (added) wasAdded = true;
			}
			
			// if it wasn't added, make a new region
			if (!wasAdded) {
				MoonShotRegion region = new MoonShotRegion(match);
				regions.add(region);
			}
			
		}
	}
	
	
	public static void printMatches(ArrayList<Match> unifiedMatches, String fileName) {
		ArrayList<String> columns = new ArrayList<String>();
		columns.add("fileLocus");
		columns.add("spectrumMD5");
		columns.add("FilePath");
		columns.add("score");
		columns.add("peptideSequence");
		columns.add("previousAminoAcid");
		columns.add("SequenceName");
		columns.add("start");
		columns.add("stop");
		columns.add("strand");
		columns.add("interest");
		columns.add("geneType");
		columns.add("geneName");
		columns.add("novelDistance");
		columns.add("novelMassDistance");
		columns.add("homologousProtein");
		columns.add("homologousPeptide");
	
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(fileName)));
			pw.println("one");
			pw.println("two");
			pw.println("three");
			//print columns
			StringBuffer line = new StringBuffer();
			for (String column: columns) {
				line.append(column + "\t");
			}
			pw.println(line);
			
			//print matches
			for (Match match: unifiedMatches) {
				if (match.getInt("novelDistance") != 0) {
					line = new StringBuffer();
					for (String column: columns) {
						line.append(match.get(column) + "\t");
					}
					pw.println(line);
				}
			}
			pw.flush();
			pw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static ArrayList<Match> getUnifiedMatches(
			File reportsFolder, 
			String novelPath, 
			boolean getGencode, 
			boolean getNovelDistance) {
		

		
		//go through each folder; load matches from our novelPath folder
		ArrayList<Match> allMatches = new ArrayList<Match>();
		for (File folder: reportsFolder.listFiles()) {
			if (folder.isDirectory()) {
				for (File novelFolder: folder.listFiles()) {
					if (novelFolder.getName().indexOf(novelPath) != -1) {
						File reportFile = new File(novelFolder, "/report.txt");
						if (reportFile.exists()) {
							ArrayList<Match> matches = Match.loadMatches(reportFile);
							allMatches.addAll(matches);
						}
					}
				}
			}
			
		}
		U.p("found: " + allMatches.size());
		
		
		
		//assign gene names to matches
		if (getGencode) {
			U.p("getting gencode regions");
			assignGencodeAnnotations(allMatches);
		}
		
		
		//foreach match, find if the novel distance is 0 or 1
		if (getNovelDistance) {
			U.p("finding novel distance");
			Sequence_Protein uniProt = new Sequence_Protein(new File("/Users/risk2/PeppyData/public/sequences/protein/UniProt-HUMAN-20130918.fasta"));
			ArrayList<Protein> proteins = uniProt.getProteinsFromDatabase(false, false);
			for(Match match: allMatches) {
				String peptide = match.getString("peptideSequence");
				int peptideMiddle = peptide.length() / 2;
				int peptideQuarter = peptide.length() / 4;
				String halfA = peptide.substring(0, peptideMiddle);
				String halfB = peptide.substring(peptideMiddle, peptide.length());
				String halfMiddle = peptide.substring(peptideQuarter, peptide.length());
				int distance = 100;
				double massDistance = -1;
				String proteinHomology = "-";
				String similarPeptide = "-";
				for(Protein protein: proteins) {
					String proteinSequence = protein.getAcidString();
					
					//if our peptide is a perfect match
					if (proteinSequence.indexOf(peptide) != -1) {
						distance = 0;
						massDistance = 0;
						proteinHomology = protein.getName();
						break;
					}
					
					//if our peptide has a strong homology
					if(
							proteinSequence.indexOf(halfA) != -1 || 
							proteinSequence.indexOf(halfB) != -1 ||
							proteinSequence.indexOf(halfMiddle) != -1) {
						int maxCommonCount = 0;
						int maxCommonCountIndex = 0;
						
						//special case for partial match near start of protein
						for (int peptideIndex = 1; peptideIndex < peptide.length(); peptideIndex++) {
							String subPeptide = peptide.substring(peptideIndex, peptide.length());
							int proteinLoc = proteinSequence.indexOf(subPeptide);
							if (proteinLoc != -1) {
								if (proteinLoc < peptide.length()) {
									maxCommonCount = peptide.length() - peptideIndex;
									maxCommonCountIndex = 0;
									break;
								}
							}
						}
						
						//slide the peptide along the protein
						char peptideAcid, proteinAcid;
						boolean acidsEqual = false;
						for (int proteinIndex = 0; proteinIndex < proteinSequence.length() - peptide.length(); proteinIndex++) {
							int commonCount = 0;
							for (int peptideIndex = 0; peptideIndex < peptide.length(); peptideIndex++) {
								peptideAcid = peptide.charAt(peptideIndex);
								proteinAcid = proteinSequence.charAt(proteinIndex + peptideIndex);
								
								//resetting acidsEqual
								acidsEqual = false;
								
								//if they are actually equal
								if (peptideAcid == proteinAcid) acidsEqual = true;
								
								//taking mass equivalence of i and l into account
								if (peptideAcid == 'I' || peptideAcid == 'L') {
									if (peptideAcid == 'I' || peptideAcid == 'L') {
										acidsEqual = true;
									}
								}
								
								//add to acids in common if equality is determined
								if (acidsEqual) {
									commonCount++;
								}
							}
							if (commonCount > maxCommonCount) {
								maxCommonCount = commonCount;
								maxCommonCountIndex = proteinIndex;
							}
						}
						
						//special case for partial match at end of protein
						for (int peptideIndex = 1; peptideIndex < peptide.length(); peptideIndex++) {
							String subPeptide = peptide.substring(0, peptideIndex);
							int endIndex = proteinSequence.length() - peptide.length();
							int proteinLoc = proteinSequence.indexOf(subPeptide, endIndex);
							if (proteinLoc != -1) {
								if (peptideIndex > maxCommonCount) {
									maxCommonCount = peptideIndex;
									maxCommonCountIndex = endIndex;
								}
							}
						}
						
						int similarPeptideStop = maxCommonCountIndex + peptide.length();
						if (similarPeptideStop > proteinSequence.length()) similarPeptideStop = proteinSequence.length();
						similarPeptide = proteinSequence.substring(maxCommonCountIndex, similarPeptideStop);
						
						double similarMass = (new Peptide(similarPeptide)).getMass();
						double originalMass = (new Peptide(peptide)).getMass();
						massDistance = originalMass - similarMass;
						proteinHomology = protein.getName();
						
						distance = peptide.length() - maxCommonCount;
						break;
					}
					
					
				}
				match.set("novelDistance", distance);
				match.set("novelMassDistance", massDistance);
				match.set("homologousProtein", proteinHomology);
				match.set("homologousPeptide", similarPeptide);
				
			}
		}
		
		return allMatches;
		
	}
	
	
	public static void remapTargetedMatches(ArrayList<Match> matches) {
		for (Match match: matches) {
			//sample protein name:
			//chr11; strand:rev; frame:2; start:27910902; stop:27910257
			String proteinName = match.getString("SequenceName");
			String [] codingComponents = proteinName.split(";");
			
			//sequence
			String sequence = codingComponents[0];
			
			//Strand: +, -
			String strand = codingComponents[1].substring(codingComponents[1].indexOf(":") + 1);
			if (strand.equals("fwd")) strand = "+";
			if (strand.equals("rev")) strand = "-";
			
			//where the protein begins
			int proteinStart = Integer.parseInt(codingComponents[3].substring(codingComponents[3].indexOf(":") + 1));
			
			//calculate the peptide coding region
			Integer start = (match.getInt("start") * 3) + proteinStart;
			Integer stop = start + (3 * match.getString("peptideSequence").length());
			
			//adjust the match for genomic coordinates
			match.set("SequenceName", ">" + sequence);
			match.set("start", start);
			match.set("stop", stop);
			match.set("strand", strand);	
		}
	}
	
	
	public static ArrayList<Region> loadGeneRegions(String regionsPath) {
		ArrayList<Region> geneRegions = new ArrayList<Region>();
		try {
			File gencodeFile = new File(regionsPath);
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
		
		return geneRegions;
	}
	
	
	public static void assignGencodeAnnotations(ArrayList<Match> allMatches) {
		ArrayList<Region> geneRegions = loadGeneRegions("resources/gencode/gencodeReduced.gtf");
		Hashtable<String, Region> foundRegions = new Hashtable<String, Region>();
		for (Match match: allMatches) {
			int startLocus = match.getInt("start");
			int stopLocus = match.getInt("stop");
			String sequence = match.getString("sequenceName");
			if (sequence.startsWith(">")) sequence = sequence.substring(1);
			for (Region region: geneRegions) {
				if (!region.getSequence().equals(sequence)) continue;
				boolean startInside = (startLocus >= region.getStart() && startLocus <= region.getStop());
				boolean stopInside = (stopLocus >= region.getStart() && stopLocus <= region.getStop());
				if (stopInside || startInside) {
					for (String tkName: Definitions.protoOncoList) {
						if (region.getName().equals(tkName)) {
							match.set("interest", "proto-oncogene"); break;
						}
					}
					for (String tkName: Definitions.gfList) {
						if (region.getName().equals(tkName)) {
							match.set("interest", "growth factor"); break;
						}
					}
					for (String tkName: Definitions.oncoList) {
						if (region.getName().equals(tkName)) {
							match.set("interest", "oncogene"); break;
						}
					}
					for (String tkName: Definitions.gtpaseActivationList) {
						if (region.getName().equals(tkName)) {
							match.set("interest", "GTPase activation"); break;
						}
					}
					for (String tkName: Definitions.tkList) {
						if (region.getName().equals(tkName)) {
							match.set("interest", "tyrosine-kinase"); break;
						}
					}
					
					match.set("geneType", region.getDescription());
					match.set("geneName", region.getName());
					
					region.addMatch(match);
					foundRegions.put(region.getName(), region);
					
					break;
				}
			}
		}
	}

}
