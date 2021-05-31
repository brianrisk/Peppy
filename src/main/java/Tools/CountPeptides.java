package Tools;

import Peppy.*;

import java.util.ArrayList;

public class CountPeptides {

    public static void main(String args[]) {
        /* set up initial state */
        Peppy.init(args);

        U.p("counting peptides for this many missed cleavages: " + Properties.numberOfMissedCleavages);
        U.p("counting the number of peptides...");

        int peptideCount = 0;

        ArrayList<Sequence> sequences = Sequence.loadSequenceFiles(Properties.sequenceDirectoryOrFile);

        /* track which sequence we are examining */
        int sequenceIndex = 0;

        /* loops until we have gone through all of our sequences */
        while (true) {


            /* Extract a decent size of peptides.  Sequences may be short, so this
             * goes through each sequence and extracts peptides until a desired
             * threshold has been been reached.
             */

            /*  Initialize our list of peptides */
            ArrayList<Peptide> peptides = new ArrayList<Peptide>(Properties.desiredPeptideDatabaseSize);
            peptides.ensureCapacity(Properties.desiredPeptideDatabaseSize);

            /* This is where we get a chunk of peptides */
            ArrayList<Peptide> peptideSegment = new ArrayList<Peptide>();

            /* if we use a region, extract that region, else go through all sequences getting a chunk at a time. */
            if (Properties.useSequenceRegion) {

                /* collect peptides */
                peptideSegment = sequences.get(sequenceIndex).extractMorePeptides(false);
                peptides.addAll(peptideSegment);

                /* reset the sequence so that the next batch of spectra can scan it */
                sequences.get(sequenceIndex).reset();

            } else {
                while (peptides.size() < Properties.desiredPeptideDatabaseSize) {

                    /* clear previous chunk of peptides and reclaim memory */
                    if (peptideSegment != null) {
                        peptideSegment.clear();
                        System.gc();
                    }

                    /* collect peptides */
                    peptideSegment = sequences.get(sequenceIndex).extractMorePeptides(false);

                    /* advance to the next sequence if we don't have any more peptides in this sequence */
                    if (peptideSegment == null) {
                        sequences.get(sequenceIndex).reset();
                        sequenceIndex++;

                        /* add peptides to the main list if we have some to add */
                    } else {
                        peptides.addAll(peptideSegment);
                    }

                    /* if we have advanced past the last sequence, then exit this loop */
                    if (sequenceIndex == sequences.size()) {
                        break;
                    }

                }
            }


            /* report */
            peptideCount += peptides.size();


            /* free up memory */
            peptides.clear();
            System.gc();

            /* break if we have covered our last sequence or if we are only using a region */
            if (sequenceIndex == sequences.size() || Properties.useSequenceRegion) {
                break;
            }

        }

        U.p("We found this many total peptides: " + peptideCount);

    }

}
