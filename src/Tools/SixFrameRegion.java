package Tools;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import Peppy.Definitions;
import Peppy.DigestionThread_DNA;
import Peppy.Peppy;
import Peppy.Properties;
import Peppy.Sequence;
import Peppy.Sequence_DNA;
import Utilities.U;

/**
 * A tool which lets you specify a region in a sequence and it comes up with
 * a FASTA file which is the six frame translation of that region
 * @author Brian Risk
 *
 */
public class SixFrameRegion {
	
	private static String sequenceDNA;
	
	public static void main(String args[]) {
		U.p("Starting...");
		
		/* initializing */
		Peppy.init(args);
		
		/* loading in the nucleotide sequence */
		ArrayList<Sequence> sequences = Sequence.loadSequenceFiles(Properties.sequenceDirectoryOrFile);
		Sequence_DNA sequenceFile = (Sequence_DNA) sequences.get(0);
		sequenceDNA = sequenceFile.getNucleotideSequences().get(0).getSequence();
		
		
		int startPosition = Properties.sequenceRegionStart;
		int stopPosition = Properties.sequenceRegionStop;
		
//		U.p(sequenceDNA.subSequence(startPosition, stopPosition));
		
		String chrName = U.getFileNameWithoutSuffix(sequenceFile.getSequenceFile());
		
		/* produce and write the frames */
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("six frame translation.txt")));
			
			/* do the forwards frames */
			for (int i = 0; i < 3; i++) {
				String frame = getFrame (startPosition + i, stopPosition, true);
				
				/* write the header */
				String header = ">" + chrName + " forwards strand: " + (startPosition + i) + ", " + stopPosition + " ";
				pw.println(header);
				
				/* print the acids, newline every 50 chars */
				for (int j = 0; j < frame.length(); j++) {
					pw.print(frame.charAt(j));
					if (j % 50 == 0 && j > 0) {
						pw.println();
					}
				}
				
				/* print extra spacer */
				pw.println();
				pw.println();
				
			}
			
			/* do the reverse frames */
			for (int i = 0; i < 3; i++) {
				String frame = getFrame (stopPosition - i, startPosition, false);
				
				/* write the header */
				String header = ">" + chrName + " reverse strand: " + (startPosition) + ", " + (stopPosition - i) + " ";
				pw.println(header);
				
				/* print the acids, newline every 50 chars */
				for (int j = 0; j < frame.length(); j++) {
					pw.print(frame.charAt(j));
					if (j % 50 == 0 && j > 0) {
						pw.println();
					}
				}
				
				/* print extra spacer */
				pw.println();
				pw.println();
				
			}
			
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		U.p("done");
		
	}
	
	private static String getFrame (int startPosition, int stopPosition, boolean isForwardsStrand) {

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
				aminoAcid = Definitions.aminoAcidList[DigestionThread_DNA.indexForCodonArray(codon, isForwardsStrand)];
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
