package Navigator;

import java.awt.Color;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.zip.GZIPOutputStream;

import Peppy.U;
import Reports.UCSC;

/**
 * Copyright 2012, Brian Risk
 * Released under the Netscape Public License
 * 
 * @author Brian Risk
 *
 */
public class SaveBedFile {
	
	String trackName;
	double scoreCutoff;
	ArrayList<Match> contaminant_protein = new ArrayList<Match>();
	ArrayList<Match> reference_protein = new ArrayList<Match>();
	ArrayList<Match> reference_genome = new ArrayList<Match>();
	ArrayList<Match> subject_genome = new ArrayList<Match>();
	ArrayList<Match> disease_genome = new ArrayList<Match>();
	
	/* keeping track of how many we find */
	private int contaminantTotal = 0;
	private int referenceProteinTotal = 0;
	private int referenceGenomeTotal = 0;
	private int subjectGenomeTotal = 0;
	private int diseaseGenomeTotal = 0;
	private int spectrumTotal = 0;
	
	final int CONTAMINANT_PROTEIN = 0;
	final int REFERENCE_PROTEIN = 1;
	final int REFERENCE_GENOME = 2;
	final int SUBJECT_GENOME = 3;
	final int DISEASE_GENOME = 4;
	
	public static void main(String arg[]) {
		U.p("loading results");
		
		SaveBedFile track;
		double threshold;
		
		/* Karen Anderson */
		threshold = 15;
		track = new SaveBedFile("2012-04", threshold);
		track.addReferenceProteinMatches("/Users/risk2/Sites/research/karen-anderson/karen-anderson 2012-04/1 2012-04 - MOUSE.fasta/report.txt");
		track.addReferenceGenomeMatches("/Users/risk2/Sites/research/karen-anderson/karen-anderson 2012-04/2 2012-04 - mouse/report.txt");
		track.save();
		
//		/* Wash U */
//		threshold = 15.09;
//		track = new SaveBedFile("WHIM 2 - 33", threshold);
//		track.addContaminantProteinMatches("/Users/risk2/PeppyData/WashU/reports/WHIM2/WashU 33 WHIM2 mouse/report.txt");
//		track.addReferenceProteinMatches("/Users/risk2/PeppyData/WashU/reports/WHIM2/WashU 33 WHIM2 human/report.txt");
//		track.addReferenceGenomeMatches("/Users/risk2/PeppyData/WashU/reports/WHIM2/WashU 33 WHIM2 hg19/report.txt");
//		track.save();
//		
//		track = new SaveBedFile("WHIM 2 - 41", threshold);
//		track.addContaminantProteinMatches("/Users/risk2/PeppyData/WashU/reports/WHIM2/WashU 41 WHIM2 mouse/report.txt");
//		track.addReferenceProteinMatches("/Users/risk2/PeppyData/WashU/reports/WHIM2/WashU 41 WHIM2 human/report.txt");
//		track.addReferenceGenomeMatches("/Users/risk2/PeppyData/WashU/reports/WHIM2/WashU 41 WHIM2 hg19/report.txt");
//		track.save();
//		
//		track = new SaveBedFile("WHIM 2 - 43", threshold);
//		track.addContaminantProteinMatches("/Users/risk2/PeppyData/WashU/reports/WHIM2/WashU 43 WHIM2 mouse/report.txt");
//		track.addReferenceProteinMatches("/Users/risk2/PeppyData/WashU/reports/WHIM2/WashU 43 WHIM2 human/report.txt");
//		track.addReferenceGenomeMatches("/Users/risk2/PeppyData/WashU/reports/WHIM2/WashU 43 WHIM2 hg19/report.txt");
//		track.save();
//		
//		track = new SaveBedFile("WHIM 16 - 33", threshold);
//		track.addContaminantProteinMatches("/Users/risk2/PeppyData/WashU/reports/WHIM16/WashU 33 WHIM16 mouse/report.txt");
//		track.addReferenceProteinMatches("/Users/risk2/PeppyData/WashU/reports/WHIM16/WashU 33 WHIM16 human/report.txt");
//		track.addReferenceGenomeMatches("/Users/risk2/PeppyData/WashU/reports/WHIM16/WashU 33 WHIM16 hg19/report.txt");
//		track.save();
//		
//		track = new SaveBedFile("WHIM 16 - 41", threshold);
//		track.addContaminantProteinMatches("/Users/risk2/PeppyData/WashU/reports/WHIM16/WashU 41 WHIM16 mouse/report.txt");
//		track.addReferenceProteinMatches("/Users/risk2/PeppyData/WashU/reports/WHIM16/WashU 41 WHIM16 human/report.txt");
//		track.addReferenceGenomeMatches("/Users/risk2/PeppyData/WashU/reports/WHIM16/WashU 41 WHIM16 hg19/report.txt");
//		track.save();
//		
//		track = new SaveBedFile("WHIM 16 - 43", threshold);
//		track.addContaminantProteinMatches("/Users/risk2/PeppyData/WashU/reports/WHIM16/WashU 43 WHIM16 mouse/report.txt");
//		track.addReferenceProteinMatches("/Users/risk2/PeppyData/WashU/reports/WHIM16/WashU 43 WHIM16 human/report.txt");
//		track.addReferenceGenomeMatches("/Users/risk2/PeppyData/WashU/reports/WHIM16/WashU 43 WHIM16 hg19/report.txt");
//		track.save();
		
//		/* Vanderbilt */
//		threshold = 13.2;
//		track = new SaveBedFile("Vanderbilt WHIM 2", threshold);
//		track.addContaminantProteinMatches("/Users/risk2/PeppyData/vanderbilt/reports/WHIM2/Vanderbilt WHIM2 mouse/report.txt");
//		track.addReferenceProteinMatches("/Users/risk2/PeppyData/vanderbilt/reports/WHIM2/Vanderbilt WHIM2 human/report.txt");
//		track.addReferenceGenomeMatches("/Users/risk2/PeppyData/vanderbilt/reports/WHIM2/Vanderbilt WHIM2 hg19/report.txt");
//		track.save();
//		
//		track = new SaveBedFile("Vanderbilt WHIM 16", threshold);
//		track.addContaminantProteinMatches("/Users/risk2/PeppyData/vanderbilt/reports/WHIM16/Vanderbilt WHIM16 mouse/report.txt");
//		track.addReferenceProteinMatches("/Users/risk2/PeppyData/vanderbilt/reports/WHIM16/Vanderbilt WHIM16 human/report.txt");
//		track.addReferenceGenomeMatches("/Users/risk2/PeppyData/vanderbilt/reports/WHIM16/Vanderbilt WHIM16 hg19/report.txt");
//		track.save();
//		
//		
//		/* JHU */
//		threshold = 23;
//		track = new SaveBedFile("JHU WHIM 2", threshold);
//		track.addContaminantProteinMatches("/Users/risk2/PeppyData/johns-hopkins/reports/WHIM2/JHU WHIM2 mouse/report.txt");
//		track.addReferenceProteinMatches("/Users/risk2/PeppyData/johns-hopkins/reports/WHIM2/JHU WHIM2 human/report.txt");
//		track.addReferenceGenomeMatches("/Users/risk2/PeppyData/johns-hopkins/reports/WHIM2/JHU WHIM2 hg19/report.txt");
//		track.save();
//		
//		track = new SaveBedFile("JHU WHIM 16", threshold);
//		track.addContaminantProteinMatches("/Users/risk2/PeppyData/johns-hopkins/reports/WHIM16/JHU WHIM16 mouse/report.txt");
//		track.addReferenceProteinMatches("/Users/risk2/PeppyData/johns-hopkins/reports/WHIM16/JHU WHIM16 human/report.txt");
//		track.addReferenceGenomeMatches("/Users/risk2/PeppyData/johns-hopkins/reports/WHIM16/JHU WHIM16 hg19/report.txt");
//		track.save();
//		
//		/* UNC */
//		threshold = 21.2;
//		track = new SaveBedFile("UNC WHIM 2", threshold);
//		track.addContaminantProteinMatches("/Users/risk2/PeppyData/UNC/reports/WHIM2/UNC WHIM2 mouse/report.txt");
//		track.addReferenceProteinMatches("/Users/risk2/PeppyData/UNC/reports/WHIM2/UNC WHIM2 human/report.txt");
//		track.addReferenceGenomeMatches("/Users/risk2/PeppyData/UNC/reports/WHIM2/UNC WHIM2 hg19/report.txt");
//		track.save();
//		
//		track = new SaveBedFile("UNC WHIM 16", threshold);
//		track.addContaminantProteinMatches("/Users/risk2/PeppyData/UNC/reports/WHIM16/UNC WHIM16 mouse/report.txt");
//		track.addReferenceProteinMatches("/Users/risk2/PeppyData/UNC/reports/WHIM16/UNC WHIM16 human/report.txt");
//		track.addReferenceGenomeMatches("/Users/risk2/PeppyData/UNC/reports/WHIM16/UNC WHIM16 hg19/report.txt");
//		track.save();
		
		
		U.p("done");
	}
	
	public SaveBedFile(String trackName, double scoreCutoff) {
		this.trackName = trackName;
		this.scoreCutoff = scoreCutoff;
	}
	
	public ArrayList<Match> loadMatches(String fileName) {
		ArrayList<Match> matches = Match.loadMatches(new File(fileName));
		int passCount = 0;
		for (Match match: matches) {
			if (match.getScore() >= scoreCutoff) passCount++;
		}
		ArrayList<Match> out = new ArrayList<Match>(passCount);
		for (Match match: matches) {
			if (match.getScore() >= scoreCutoff) out.add(match);
		}
		return out;
	}
	
	public void addContaminantProteinMatches(String fileName) {
		ArrayList<Match> matches = loadMatches(fileName);
		for (Match match: matches) match.set("label", CONTAMINANT_PROTEIN);
		for (Match match: matches) match.set("isProtein", true);
		contaminant_protein.addAll(matches);
	}
	
	public void addReferenceProteinMatches(String fileName) {
		ArrayList<Match> matches = loadMatches(fileName);
		for (Match match: matches) match.set("label", REFERENCE_PROTEIN);
		for (Match match: matches) match.set("isProtein", true);
		reference_protein.addAll(matches);
	}
	
	public void addReferenceGenomeMatches(String fileName) {
		ArrayList<Match> matches = loadMatches(fileName);
		for (Match match: matches) match.set("label", REFERENCE_GENOME);
		for (Match match: matches) match.set("isProtein", false);
		reference_genome.addAll(matches);
	}
	
	public void addSubjectGenomeMatches(String fileName) {
		ArrayList<Match> matches = loadMatches(fileName);
		for (Match match: matches) match.set("label", SUBJECT_GENOME);
		for (Match match: matches) match.set("isProtein", false);
		subject_genome.addAll(matches);
	}
	
	public void addDiseaseGenomeMatches(String fileName) {
		ArrayList<Match> matches = loadMatches(fileName);
		for (Match match: matches) match.set("label", DISEASE_GENOME);
		for (Match match: matches) match.set("isProtein", false);
		disease_genome.addAll(matches);
	}
	
	
	public void save() {
		
		File folder = new File("BED files");
		folder.mkdir();
		
		/* the best matches */
		MatchTable best = new MatchTable(true);
		for (Match match: contaminant_protein) best.put(match.getString("spectrumMD5"), match);
		for (Match match: reference_protein) best.put(match.getString("spectrumMD5"), match);
		for (Match match: reference_genome) best.put(match.getString("spectrumMD5"), match);
		for (Match match: subject_genome) best.put(match.getString("spectrumMD5"), match);
		for (Match match: disease_genome) best.put(match.getString("spectrumMD5"), match);
		
		U.p("best has this many elements: " + best.getHashtable().size());

		
		/* save bed file matches */
		U.p("saving bedfile");
		try {
			File trackFile = new File(folder, trackName + " bed.txt");
			PrintWriter bedWriter = new PrintWriter(new BufferedWriter(new FileWriter(trackFile)));
			
			PrintWriter mouseWriter = new PrintWriter(new BufferedWriter(new FileWriter(new File(folder, trackName + "mouseOnly.html"))));
			mouseWriter.print("<html><body><ul>");
			
			PrintWriter tumorWriter = new PrintWriter(new BufferedWriter(new FileWriter(new File(folder, trackName +"tumorOnly.html"))));
			tumorWriter.print("<html><body><ul>");
			
			
			
			/* define the track properties */
			bedWriter.println("track name=\"" + trackName +"\" description=\"" + trackName + "\" visibility=2 itemRgb=\"On\"");
			
			/* setting up variables */
			Hashtable<String, ArrayList<Match>> bestHashtable = best.getHashtable();
			
			int label;
			String md5;
			ArrayList<Match> spectrumMatches ;
			
			
			Enumeration<String> e = bestHashtable.keys();
			while (e.hasMoreElements()) {
				md5 = e.nextElement();

				spectrumMatches = bestHashtable.get(md5);
				
				/* these variables will mark where we are finding these matches.
				 * they will be used further down 
				 */
				boolean isInProtein = false;
				boolean isInReferenceProtein = false;
				boolean isInReferenceGenome = false;
				boolean isInSubjectGenome = false;
				
				/*
				 * Finding out if the best match can be found in a protein, and
				 * which specific databases contain it.
				 */
				for (Match match: spectrumMatches) {
					/* is it a protein */
					if (match.getBoolean("isProtein")) {
						isInProtein = true;
					}
					
					/* what is the label? */
					label = match.getInt("label");
					if (label != CONTAMINANT_PROTEIN) {
						if (label == REFERENCE_PROTEIN) {
							isInReferenceProtein = true;
						} else {
							if (label == REFERENCE_GENOME) {
								isInReferenceGenome = true;
							} else {
								if (label == SUBJECT_GENOME) {
									isInSubjectGenome = true;
								}
							}
						}
					}
				}
				
				/* the dominant label is how these tracks will be colored. */
				spectrumTotal++;
				int dominantLabel = REFERENCE_PROTEIN;
				if (isInProtein) {
					if (!isInReferenceProtein) {
						dominantLabel = CONTAMINANT_PROTEIN;
						contaminantTotal++;
					} else {
						referenceProteinTotal++;
					}
				} else {
					if (isInReferenceGenome) {
						dominantLabel = REFERENCE_GENOME;
						referenceGenomeTotal++;
					} else {
						if (isInSubjectGenome) {
							dominantLabel = SUBJECT_GENOME;
							subjectGenomeTotal++;
						} else {
							dominantLabel = DISEASE_GENOME;
							diseaseGenomeTotal++;
						}
					}
				}
				
				
				
				/* delete duplicates.  If a peptide appears in a reference protein and genome and
				 * subject genome, we don't want that same peptide appearing 3 times
				 * 
				 * Also removes the protein matches as they do not have genome coordinates
				 * 
				 * Mildly confusing.  Sorry about this.
				 * 
				 * Even if the dominant label is from the reference or contaminant protein,
				 * we are keeping just the reference genome hits.  Again, this is because these
				 * are the matches that have genomic coordinates.
				 * 
				 * Then, we make special exceptions for the subject and disease genomes
				 */
				ArrayList<Match> reducedSpectrumMatches = new ArrayList<Match>();
				int keepLabel = REFERENCE_GENOME;
				if (dominantLabel == SUBJECT_GENOME)  keepLabel = SUBJECT_GENOME;
				if (dominantLabel == DISEASE_GENOME)  keepLabel = DISEASE_GENOME;
				for (Match match: spectrumMatches) {
					if (match.getInt("label") == keepLabel) reducedSpectrumMatches.add(match);
				}
				spectrumMatches = reducedSpectrumMatches;
				
				
				/* printing our BED reports */
				for (Match match: spectrumMatches) {
					String link = UCSC.getLink(match.getInt("start"), match.getInt("stop"), match.getString("SequenceName"));
					if (dominantLabel == CONTAMINANT_PROTEIN) {
						bedWriter.println(getBedLine(match, Color.orange));
						mouseWriter.print("<li><a href=\"" + link + "\">" + match.getString("peptideSequence") +"</a></li>" );
					} else {
						if (dominantLabel == REFERENCE_PROTEIN) {
							bedWriter.println(getBedLine(match, Color.BLACK));
						} else {
							if (dominantLabel == REFERENCE_GENOME) {
								bedWriter.println(getBedLine(match, Color.BLUE));
							} else {
								if (dominantLabel == SUBJECT_GENOME) {
									bedWriter.println(getBedLine(match, Color.RED));
								} else {
									bedWriter.println(getBedLine(match, Color.GREEN));
									tumorWriter.print("<li><a href=\"" + link + "\">" + match.getString("peptideSequence") +"</a></li>" );
								}
							}
						} 
							
					}
					
				}
				
				
			}
			
			
			
			bedWriter.flush();
			bedWriter.close();
			
			mouseWriter.print("</ul></body></html>");
			mouseWriter.flush();
			mouseWriter.close();
			
			
			tumorWriter.print("</ul></body></html>");
			tumorWriter.flush();
			tumorWriter.close();
			
			

			/* gzip the file */
			String gzipFileName = (trackName + " bed.gz").replace(' ', '-');
			GZIPOutputStream gzip = new GZIPOutputStream(new BufferedOutputStream(new FileOutputStream(new File(folder, gzipFileName))));
			int nugget;
			BufferedInputStream bis = new BufferedInputStream(new FileInputStream(trackFile));
			nugget = bis.read();
			while (nugget != -1) {
				gzip.write(nugget);
				nugget = bis.read();
			}
			bis.close();
			
			gzip.finish();
			gzip.flush();
			gzip.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static String getBedLine(Match match, Color color) {
		color = getShade(match, color);
		int thickStart = match.getInt("start");
		if (match.getBoolean("isModified")) {
			thickStart = (thickStart + match.getInt("stop")) / 2;
		}
		StringBuffer line = new StringBuffer();
		line.append(match.getString("SequenceName"));
		line.append(" ");
		line.append(match.getInt("start"));
		line.append(" ");
		line.append(match.getInt("stop"));
		line.append(" ");
		line.append(match.getString("peptideSequence"));
		line.append(" ");
		line.append("0");
		line.append(" ");
		line.append(match.getString("Strand"));
		line.append(" ");
		line.append(thickStart);
		line.append(" ");
		line.append(match.getInt("stop"));
		line.append(" ");
		line.append(color.getRed());
		line.append(",");
		line.append(color.getGreen());
		line.append(",");
		line.append(color.getBlue());
		return line.toString();
	}
	
	private static Color getShade(Match match, Color color) {

	
		int shade = 0;
		double score = match.getScore();
		
//		if (score > 23.97065894069367) shade = 5;
//		if (score > 25.620117028827554) shade = 4;
//		if (score > 26.095221891211644) shade = 3;
//		if (score > 27.380266416210844) shade = 2;
//		if (score > 27.614889698140427) shade = 1;
//		if (score > 29.00520687654409) shade = 0;
		
		int red = colorNumberCorrect((shade * 255 + color.getRed()) / (shade + 1));
		int green = colorNumberCorrect((shade * 255 + color.getGreen()) / (shade + 1));
		int blue = colorNumberCorrect((shade * 255 + color.getBlue()) / (shade + 1));
		return new Color(red, green, blue);
	}
	
	private static int colorNumberCorrect(int number) {
		if (number < 0) return 0;
		if (number > 255) return 255;
		return number;
	}
	
	
	/**
	 * returns setA minus setB
	 * @param setA
	 * @param setB
	 * @return
	 */
	public static ArrayList<Match> minus(ArrayList<Match> setA, ArrayList<Match> setB) {
		ArrayList<Match> out = new ArrayList<Match>();
		for (Match setAMatch: setA) {
			boolean found = false;
			for (Match setBMatch: setB) {
				if (setAMatch.equals(setBMatch)) {
					found = true;
					break;
				}
			}
			if (! found) {
				out.add(setAMatch);
			}
		}
		return out;
	}
	


}
