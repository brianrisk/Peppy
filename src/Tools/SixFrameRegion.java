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
	
	private static int minimumProteinLength = 8;
	
	
	public static void main(String args[]) {
		U.p("Starting six frame region translation...");
		
		/* initializing */
		Peppy.init(args);
		
		/* sequence will always be DNA */
		Properties.isSequenceFileDNA = true;
 		
//		createDatabase(Properties.sequenceDirectoryOrFile, Properties.sequenceRegionStart, Properties.sequenceRegionStop);
		
//		File interestLocations = new File("regions-of-interest.csv");
//		createDatabasesFromFile(interestLocations, "/Users/risk2/PeppyData/public/sequences/dna/HG19/");
		
//		createDatabaseOfGenome("/Users/risk2/PeppyData/public/sequences/dna/HG19/", "hg19");
//		createDatabaseOfGenome("/Users/risk2/PeppyData/WashU/sequences/whim16/dna/germline/", "w16-germline");
//		createDatabaseOfGenome("/Users/risk2/PeppyData/WashU/sequences/whim16/dna/tumor/", "w16-tumor");
//		createDatabaseOfGenome("/Users/risk2/PeppyData/WashU/sequences/whim16/dna/xeno/", "w16-xeno");
//		createDatabaseOfGenome("/Users/risk2/PeppyData/WashU/sequences/whim2/dna/germline/", "w2-germline");
//		createDatabaseOfGenome("/Users/risk2/PeppyData/WashU/sequences/whim2/dna/tumor/", "w2-tumor");
//		createDatabaseOfGenome("/Users/risk2/PeppyData/WashU/sequences/whim2/dna/xeno/", "w2-xeno");
		
//		createDatabase(new File("/Users/risk2/PeppyData/public/sequences/dna/HG19/chr17.fa"), 26671724, 26675221, "seleno-regions");
//		createDatabase(new File("/Users/risk2/PeppyData/public/sequences/dna/HG19/chr3.fa"), 142657639, 142671525, "seleno-regions");
//		createDatabase(new File("/Users/risk2/PeppyData/public/sequences/dna/HG19/chr2.fa"), 178495137, 178528740, "seleno-regions");
//		
		
		createDatabase(new File("/Users/risk2/PeppyData/public/sequences/dna/HG19/chr7.fa"), 48000000, 49000000, "GM-chr7-region");
		
		
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
		int windowRadius = 50000;
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(interestLocations));
			String line = br.readLine();
			while (line != null) {
				String [] chunks = line.split(",");
				String acidSequence = chunks[0];
				String chrName = chunks[1];
				int locus = Integer.parseInt(chunks[2]);
				File sequenceFile = new File(sequenceDirectory, chrName + ".fa");
				int startPosition = locus - windowRadius;
				int stopPosition = locus + windowRadius;
				createDatabase(sequenceFile, startPosition, stopPosition, destinationFolderName);
				/* read the next line */
				line = br.readLine();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	
	public static void createDatabaseOfGenome(String sequenceDirectoryString, String destinationFolderName) {
		/* let the user see which one we're on */
		U.p(destinationFolderName);
		
		ArrayList<Sequence> sequences = Sequence.loadSequenceFiles(sequenceDirectoryString);	
		for (Sequence sequence: sequences) {
			U.p(sequence.getSequenceFile().getName());
			createDatabase(sequence.getSequenceFile(), -1, -1, destinationFolderName);
		}
	}
	
		
	/**
	 * setting startPosition to -1 will make it do the whole sequence
	 * @param sequence
	 * @param startPosition
	 * @param stopPosition
	 */
	public static void createDatabase(File sequence, int startPosition, int stopPosition, String destinationFolderName) {
		
		boolean wholeSequence = false;
		if (startPosition == -1) {
			wholeSequence = true;
		}
		
		ArrayList<Sequence> sequences = Sequence.loadSequenceFiles(sequence);
		Sequence_DNA sequenceFile = (Sequence_DNA) sequences.get(0);
		
		String chrName = U.getFileNameWithoutSuffix(sequenceFile.getSequenceFile());
		
		/* where we store the fasta header info */
		String proteinHeader;
		
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
			
			/* do the forwards frames */
			for (int i = 0; i < 3; i++) {
				String frame = getFrame (startPosition + i, stopPosition, true, sequenceDNA);
				
				/* string buffer and string to hold the protein */
				StringBuffer sb = new StringBuffer();
				String protein;
				
				/* where we hold the position */
				int proteinStart, proteinStop;
				
				/* print the proteins */
				for (int j = 0; j < frame.length(); j++) {
					if (frame.charAt(j) == '.') {
						
						protein = sb.toString();
						
						/* print if protein */
						if (protein.length() > minimumProteinLength) {
							proteinStart =  (startPosition + i + (j * 3));
							proteinStop = (startPosition + i + (j * 3) + (protein.length() * 3));
							
							proteinHeader = ">" + chrName + "; strand:fwd; frame:" + i +"; start:" + proteinStart + "; stop:" + proteinStop;
							pw.println(proteinHeader );
							StringBuffer lineBuffer = new StringBuffer();
							for (int aaIndex = 0; aaIndex < protein.length(); aaIndex++) {
								lineBuffer.append(protein.charAt(aaIndex));
								if (aaIndex % 80 == 79) {
									pw.println(lineBuffer.toString());
									lineBuffer = new StringBuffer();
								}
							}
							if (lineBuffer.length() > 0) {
								pw.println(lineBuffer.toString());
							}
							pw.println();
						}
						
						/* clear out the string buffer */
						sb = new StringBuffer();
						
					} else {
						sb.append(frame.charAt(j));
					}
				}
				
				
			}
			
			/* do the reverse frames */
			for (int i = 0; i < 3; i++) {
				String frame = getFrame (stopPosition - i, startPosition, false, sequenceDNA);
				
				/* string buffer and string to hold the protein */
				StringBuffer sb = new StringBuffer();
				String protein;
				
				/* where we hold the position */
				int proteinStart, proteinStop;
				
				/* print the proteins */
				for (int j = 0; j < frame.length(); j++) {
					if (frame.charAt(j) == '.') {
						
						protein = sb.toString();
						
						/* print if protein */
						if (protein.length() > minimumProteinLength) {
							proteinStart =  stopPosition - ( i + (j * 3));
							proteinStop = stopPosition - ( i + (j * 3) + (protein.length() * 3));
							proteinHeader = ">" + chrName + "; strand:rev; frame:" + i +"; start:" + proteinStart + "; stop:" + proteinStop;
							pw.println(proteinHeader);
							StringBuffer lineBuffer = new StringBuffer();
							for (int aaIndex = 0; aaIndex < protein.length(); aaIndex++) {
								lineBuffer.append(protein.charAt(aaIndex));
								if (aaIndex % 80 == 79) {
									pw.println(lineBuffer.toString());
									lineBuffer = new StringBuffer();
								}
							}
							if (lineBuffer.length() > 0) {
								pw.println(lineBuffer.toString());
							}
							pw.println();
						}
						
						/* clear out the string buffer */
						sb = new StringBuffer();
						
					} else {
						sb.append(frame.charAt(j));
					}
				}
				
			}
			
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	
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
