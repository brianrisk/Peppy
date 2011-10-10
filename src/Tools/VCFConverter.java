package Tools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import Peppy.Nucleotides;
import Peppy.Properties;
import Peppy.Sequence;
import Peppy.Sequence_DNA;
import Utilities.U;

public class VCFConverter {

	public static void main (String[] args) {
		U.p("Starting...");
		Peppy.Peppy.init(args);
		ArrayList<Sequence> sequences = Sequence.loadSequenceFiles(Properties.sequenceDirectoryOrFile);
		U.p("we are working with " + sequences.size() + " sequence files");
				
		Sequence_DNA sequenceFile;
		Nucleotides nucleotides;
		String nucleotideSequenceString;
		char[] bases;
		
		VCFFile vcf = new VCFFile(Properties.VCFFileString);
		
		ArrayList<VCFEntry> entries = vcf.getEntries();
		
		String chromosomeNumberString;
		
		File parentDir = new File(System.getenv("HOME")+"/peppy/mutatedGenome");
		parentDir.mkdirs();
		
		for(Sequence s : sequences) {
			sequenceFile = (Sequence_DNA) s;
			nucleotides = sequenceFile.getNucleotideSequences().get(0);
			nucleotideSequenceString = nucleotides.getSequence();
			bases = nucleotideSequenceString.toCharArray();
			chromosomeNumberString = nucleotides.getChromosomeNumberString();
			U.p("Working on sequence file " + s.getSequenceFile().getName() );
			for (VCFEntry entry : entries) {
				int startIndex = entry.getPos()-1;
				if (entry.getChrom().equals(chromosomeNumberString)) {
					for (int i = 0; i < entry.getAlt().length(); i++) {
						if (bases[startIndex+i] != entry.getRef().charAt(i)) {
							U.p("What the hell?");
						}
						bases[startIndex+i] = entry.getAlt().charAt(i);
					}
				}
			}
			
			PrintWriter pw;
			try {
				pw = new PrintWriter(new BufferedWriter (new FileWriter(new File(parentDir, "chr"+chromosomeNumberString+".fa"))));
				pw.println(nucleotides.getSequenceDescription());
				for (int i = 0; i < bases.length; i++) {
					pw.print(bases[i]);
					if (i % 50 == 0 && i > 0) {
						pw.println();
					}
				}
				pw.flush();
				pw.close();
			} catch (IOException e) {
				e.printStackTrace();
			}	
		}
		U.p("... done.");
	}
}
