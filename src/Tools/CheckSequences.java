package Tools;

import java.util.ArrayList;

import Peppy.Peptide;
import Peppy.Properties;
import Peppy.Sequence;
import Utilities.U;

/**
 * This is a little tool to see if a certain peptide is appearing in the digested database
 * @author Brian Risk
 *
 */
public class CheckSequences {
	
	public static void checkSequences() {
		U.p("making sure we have the correct start and stop locations");
		Properties.isSequenceFileDNA = true;
		Properties.numberOfMissedCleavages = 1;
		ArrayList<Sequence> sequences = Sequence.loadSequenceFiles(Properties.sequenceDirectoryOrFile);
		ArrayList<Peptide> peptides;
		
		for (Sequence sequence: sequences) {
			U.p("working on sequence: " +sequence.getSequenceFile().getName());
			
			peptides = sequence.extractMorePeptides(false);
			//continually extract peptides from the sequence until there aren't anymore
			while (peptides != null) {
				for (Peptide peptide: peptides) {
					if (peptide.equals("IQDKEGIPPDQQR")) {
						U.p (peptide.getAcidSequenceString() + ": " + peptide.getStartIndex() + "\t" + peptide.getStopIndex() + "\t" + sequence.getSequenceFile().getName());
					}
				}
				peptides.clear();
				System.gc();
				peptides = sequence.extractMorePeptides(false);
			}
		}
		U.p("done");
		
	}

}
