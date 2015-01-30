package Tools.VCFConverter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import Peppy.NucleotideSequence;
import Peppy.Properties;
import Peppy.Sequence;
import Peppy.SequenceNucleotide;
import Peppy.U;

public class VCFGenomeAssembler {

	public static void main (String[] args) {
		U.p("Starting...");
		Peppy.Peppy.init(args);
		ArrayList<Sequence> sequences = Sequence.loadSequenceFiles(Properties.sequenceDirectoryOrFile);
		U.p("we are working with " + sequences.size() + " sequence files");
				
		SequenceNucleotide sequenceFile;
		NucleotideSequence nucleotideSequence;
		String nucleotideSequenceString;
		char[] bases;
		
		
		String chromosomeNumberString;
		
		File parentDir = new File("altered genome");
		parentDir.mkdirs();
		
		/* loading the VCF contents */
		VCFFile vcf = new VCFFile(Properties.VCFFileString);
		ArrayList<VCFEntry> entries = vcf.getEntries();
		
		/* report on how the loading went */
		U.p("found this many entries: " + entries.size());
		U.p("found this many invalid entries: "+ vcf.getInvalidEntries().size());
		
		/* exit if nothing */
		if (entries.size() == 0) {
			U.p("nothing to do");
			return;
		}
		
		/* track how many times the expected nucleotide was wrong */
		int numberOfTimesExpectedNucleotideWasWrong = 0;
		int numberOfTimesExpectedEqualsNew = 0;
		
		for(Sequence s : sequences) {
			sequenceFile = (SequenceNucleotide) s;
			nucleotideSequence = sequenceFile.getNucleotideSequences().get(0);
			nucleotideSequenceString = nucleotideSequence.getSequence();
			bases = nucleotideSequenceString.toCharArray();
			chromosomeNumberString = nucleotideSequence.getChromosomeNumberString();
			U.p("Working on sequence file " + s.getSequenceFile().getName() );
			for (VCFEntry entry : entries) {
				int startIndex = entry.getPos()-1;
				if (entry.getChrom().equals(chromosomeNumberString)) {
					for (int i = 0; i < entry.getAlt().length(); i++) {
						if (bases[startIndex+i] != entry.getRef().charAt(i)) {
//							U.p("What the hell?");
							numberOfTimesExpectedNucleotideWasWrong++;
						}
						if (bases[startIndex+i] == entry.getAlt().charAt(i)) {
							numberOfTimesExpectedEqualsNew++;
						}
						
						bases[startIndex+i] = entry.getAlt().charAt(i);
					}
				}
			}
			
			PrintWriter pw;
			try {
				pw = new PrintWriter(new BufferedWriter (new FileWriter(new File(parentDir, "chr"+chromosomeNumberString+".fa"))));
				pw.println(nucleotideSequence.getSequenceDescription());
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
		double percentWrong = numberOfTimesExpectedNucleotideWasWrong / entries.size();
		U.p("the expected nucleotide was wrong this many times: " + numberOfTimesExpectedNucleotideWasWrong + "(" + Properties.percentFormat.format(percentWrong) + ")");
		U.p("the expected nucleotide equaled the new one: " + numberOfTimesExpectedEqualsNew);
		U.p("... done.");
	}
}
