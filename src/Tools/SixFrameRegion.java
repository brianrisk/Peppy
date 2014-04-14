package Tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Random;

import Navigator.Region;
import Peppy.AminoAcids;
import Peppy.DigestionThread_DNA;
import Peppy.Peppy;
import Peppy.Properties;
import Peppy.Sequence;
import Peppy.Sequence_DNA;
import Peppy.U;

/**
 * A tool which lets you specify a region in a sequence and it comes up with
 * a FASTA file which is the six frame translation of that region
 * @author Brian Risk
 *
 */
public class SixFrameRegion {
	
	private static int minimumProteinLength = 0;
	private static String [] chrSuffixes = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"};
	
	public static void main(String args[]) {
		U.p("Starting six frame region translation...");
		
		/* initializing */
		Peppy.init(args);
		
		/* sequence will always be DNA */
		Properties.isSequenceFileDNA = true;

		
//		createDatabaseFromGTF();
		
		createDatabasesFromFile(
				new File("/Users/risk2/Documents/workspace/JavaGFS/targeted regions file.txt"),
				"/Users/risk2/PeppyData/public/sequences/dna/HG19",
				"CompRef regions");
		
		
		//http://www.bioline.com/calculator/01_13.html
//		String dna = "";
//		Random random = new Random();
//		char [] nucleotides = {'A', 'T', 'G', 'C'};
//		for (int i = 0; i < 1000; i++) {
//			dna += nucleotides[random.nextInt(4)];
//		}
//		U.p(dna);
//		U.p();
//		U.p(getProteins("demo", 0, dna.length(), dna));
		
		
		U.p("done");
	}
	
	
	/**
	 * This takes the chromosome and start position and creates databases with 100,000 nt diameter windows
	 * 
	 * file format:  CSV
	 * column 1: peptide sequence
	 * column 2: chromosome
	 * column 3: location
	 * 
	 * assumes:
	 * suffix is ".fa"
	 * 
	 * @param interestLocations
	 */
	public static void createDatabasesFromFile(File interestLocations, String sequenceDirectoryString, String destinationFolderName) {
		
		File sequenceDirectory = new File(sequenceDirectoryString);
		int windowRadius = 10000;
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(interestLocations));
			String line = br.readLine();
			while (line != null) {
				String [] chunks = line.split("\t");
				String acidSequence = chunks[0];
				String chrName = chunks[1];
				int locus = Integer.parseInt(chunks[2]);
				File sequenceFile = new File(sequenceDirectory, chrName + ".fa");
				int startPosition = locus - windowRadius;
				int stopPosition = locus + windowRadius;
				U.p(acidSequence);
				createDatabase(sequenceFile, startPosition, stopPosition, destinationFolderName);
				/* read the next line */
				line = br.readLine();
			}
			
			br.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public static void createDatabaseFromGTF() {
		
		/* load GTF */
		ArrayList<Region> regions = new ArrayList<Region>();
		U.p("loading gencode transcripts");
		try {
//			File gencodeFile = new File("resources/gencode/gencodeUTRReduced.gtf");
			File gencodeFile = new File("gencodeUTRReduced.gtf");
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
				regions.add(region);
				line = gencodeReader.readLine();
			}
			gencodeReader.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		/* add padding to each region */
		for (Region region: regions) {
			region.addPadding(60);
		}
		
		/* reduce regions to non-intersecting */
		ArrayList<Region> nonIntersectingRegions = new ArrayList<Region>();
		for (Region region: regions) {
			boolean intersected = false;
			for (Region regionB: nonIntersectingRegions) {
				intersected = regionB.addRegion(region);
				if (intersected) break;
			}
			if (!intersected) nonIntersectingRegions.add(region);
		}
		
		/* tally total coverage and report it */
		int coverage = 0;
		for (Region region: nonIntersectingRegions) {
			coverage += region.getCoverage();
		}
		U.p("total nucleotide coverage: " + coverage);
		U.p("original region count: " + regions.size());
		U.p("reduced cound: " + nonIntersectingRegions.size());
		//560,429,705
		
		/* produce six-frame translation of regions */
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("GTFSixFrame.fasta")));
			
			//loop through each of the chromosomes
			for (String chrSuffix: chrSuffixes) {
				
				//our chromosome name
				String chrName = "chr" + chrSuffix;
				U.p(chrName);
				
				
				//load DNA string for this chromosome
				File chrFile = new File("/Users/risk2/PeppyData/public/sequences/dna/HG19/" + chrName + ".fa");
				ArrayList<Sequence> sequences = Sequence.loadSequenceFiles(chrFile);
				Sequence_DNA sequenceFile = (Sequence_DNA) sequences.get(0);
				String sequenceDNA = sequenceFile.getNucleotideSequences().get(0).getSequence();
				
				//create translations of all regions in this chromosome
				for(Region region: nonIntersectingRegions) {
					if (region.getSequence().equals(chrName)) {
						pw.println(getProteins(chrName, region.getStart(), region.getStop(), sequenceDNA));
					}
				}
			}
			
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public static void createDatabase(File sequence, int startPosition, int stopPosition, String destinationFolderName) {
		boolean wholeSequence = (startPosition == -1);
		
		ArrayList<Sequence> sequences = Sequence.loadSequenceFiles(sequence);
		Sequence_DNA sequenceFile = (Sequence_DNA) sequences.get(0);
		
		String chrName = U.getFileNameWithoutSuffix(sequenceFile.getSequenceFile());
		
		/* the actual nucleotide sequence */
		String sequenceDNA = sequenceFile.getNucleotideSequences().get(0).getSequence();
		
		/* if whole sequence, reset start and stop now that we have the sequence loaded */
		if (wholeSequence) {
			startPosition = 0;
			stopPosition = sequenceDNA.length() - 1;
		}
		
		/* produce and write the frames */
		try {
			File destinationFolder = new File("sixFrameSequences/" + destinationFolderName);
			destinationFolder.mkdirs();
			String fileName;
			if (wholeSequence) {
				fileName = chrName +".fasta";
			} else {
				fileName = chrName + " " + startPosition + "-" + stopPosition +".fasta";
			}
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(destinationFolder, fileName))));
			
			pw.println(getProteins(chrName, startPosition, stopPosition, sequenceDNA));
			
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	private static String getProteins(String chrName, int startPosition, int stopPosition, String sequenceDNA) {
		StringBuffer out = new StringBuffer();
		
		/* do the forwards, then reverse strands */
		for(int strandIndicator = 0; strandIndicator <= 1; strandIndicator++) {
			
			//specifying which strand
			boolean isForwards = true;
			if (strandIndicator == 1) isForwards = false;
			
			//loop through our frames for this strand
			for (int frameNumber = 0; frameNumber < 3; frameNumber++) {
				String translatedFrame = "";
				
				if (isForwards) {
					translatedFrame = getFrame (startPosition + frameNumber, stopPosition, isForwards, sequenceDNA);
				} else {
					translatedFrame = getFrame (stopPosition -1 - frameNumber, startPosition - 1, isForwards, sequenceDNA);
				}
				
				/* string buffer and string to hold the protein */
				StringBuffer sb = new StringBuffer();
				
				/* print the proteins */
				for (int translatedFramePosition = 0; translatedFramePosition < translatedFrame.length(); translatedFramePosition++) {
					if (translatedFrame.charAt(translatedFramePosition) == '.') {
						
						/* append the formatted protein */
						out.append(getProtein(sb, startPosition, stopPosition, frameNumber, translatedFramePosition, chrName, isForwards));
						
						/* clear out the string buffer */
						sb = new StringBuffer();
						
					} else {
						sb.append(translatedFrame.charAt(translatedFramePosition));
					}
				}
				
				/* append the final protein protein (as sequence probably did no end with a ".") */
				out.append(getProtein(sb, startPosition, stopPosition, frameNumber, translatedFrame.length(), chrName, isForwards));
				
			}
		}
		
		return out.toString();
	}
	
	private static StringBuffer getProtein (
			StringBuffer sb, 
			int startPosition, 
			int stopPosition, 
			int frameNumber, 
			int translatedFramePosition, 
			String chrName,
			boolean isForwardsStrand) {
		
		StringBuffer out = new StringBuffer();
		
		String strand = "+";
		if (!isForwardsStrand) strand = "-";
		
		/* print if protein */
		if (sb.length() > minimumProteinLength) {
			
			/* where we hold the boundary locations */
			int proteinStop = 0, proteinStart = 0;
			
			/* determining the protein boundary locations */
			if (isForwardsStrand) {
				proteinStart = startPosition + ( frameNumber + (translatedFramePosition * 3) - (sb.length() * 3));
				proteinStop  = startPosition + ( frameNumber + (translatedFramePosition * 3));
			} else {
				proteinStop  =  stopPosition - ( frameNumber + (translatedFramePosition * 3) - (sb.length() * 3));
				proteinStart =  stopPosition - ( frameNumber + (translatedFramePosition * 3));
			}
			out.append(">" + chrName + "; strand:" + strand + "; frame:" + frameNumber +"; start:" + proteinStart + "; stop:" + proteinStop + "\r");
			StringBuffer lineBuffer = new StringBuffer();
			for (int aaIndex = 0; aaIndex < sb.length(); aaIndex++) {
				lineBuffer.append(sb.charAt(aaIndex));
				if (aaIndex % 80 == 79) {
					lineBuffer.append("\r");
					out.append(lineBuffer.toString());
					lineBuffer = new StringBuffer();
				}
			}
			if (lineBuffer.length() > 0) {
				lineBuffer.append("\r");
				out.append(lineBuffer.toString());
			}
			out.append("\r\r");
		}
		
		return out;
	}
	
	/**
	 * Does a straight translation of a DNA string
	 * 
	 * @param startPosition
	 * @param stopPosition
	 * @param isForwardsStrand
	 * @param sequenceDNA
	 * @return
	 */
	private static String getFrame (int startPosition, int stopPosition, boolean isForwardsStrand, String sequenceDNA) {

		/* translating */
		char [] codon = new char[3];
		char aminoAcid;
		int mod = 0;
		StringBuffer buildingProtein = new StringBuffer();

		int increment = 1;
		if (!isForwardsStrand) increment = -1;
		int index;
		for (index = startPosition; index != stopPosition; index += increment) {
			codon[mod] = sequenceDNA.charAt(index);
			if (mod == 2) {
				aminoAcid = AminoAcids.aminoAcidList[DigestionThread_DNA.indexForCodonArray(codon, isForwardsStrand)];
				buildingProtein.append(aminoAcid);
				
				/* reset mod */
				mod = 0;
			} else {
				mod++;
			}
		}
		
		return buildingProtein.toString();
	}

}
