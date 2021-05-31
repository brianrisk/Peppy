package Experimental;

import Peppy.*;

import java.io.File;
import java.util.ArrayList;
import java.util.Hashtable;

/**
 * Methods used to get data for the peppy score paper
 *
 * @author Brian Risk
 */
public class PaperExperiments {

    static String propertyFileString = null;
    static ArrayList<Peptide> peptides;

    public static boolean p1 = true;
    public static boolean p2 = true;
    public static boolean p3 = true;

    public static boolean p2a = true;
    public static boolean p2b = true;
    public static boolean p2c = true;
    public static boolean p2d = true;

    static double fragmentTolerance = 10;


    public static void main(String args[]) {
        /* set up initial state */
        Peppy.init(args);

//		Sequence_Protein proteinFile = new Sequence_Protein(new File("/Users/risk2/PeppyData/public/sequences/protein/UniProt_Human_2012_03.fasta"));
//		peptides = proteinFile.extractAllPeptides(false);

        int increment = 2;

        /* do the work */
        U.p("p1, p2, p3");
        p1 = true;
        p2 = true;
        p3 = true;
        for (int i = 10; i < 151; i += increment) {
            fragmentTolerance = i;
            runJobs(args);
        }

        U.p();
        U.p("p1, p3");
        p1 = true;
        p2 = false;
        p3 = true;
        for (int i = 10; i < 151; i += increment) {
            fragmentTolerance = i;
            runJobs(args);
        }

        U.p();
        U.p("p1, p2");
        p1 = true;
        p2 = true;
        p3 = false;
        for (int i = 10; i < 151; i += increment) {
            fragmentTolerance = i;
            runJobs(args);
        }

        U.p();
        U.p("p1");
        p1 = true;
        p2 = false;
        p3 = false;
        for (int i = 10; i < 151; i += increment) {
            fragmentTolerance = i;
            runJobs(args);
        }


    }


    /**
     * There may be a multiplicity of sequence and spectral directories.
     * This iterates through them in every combination.
     *
     * @param args
     */
    public static void runFunnelPeppy(String[] args) {

        Properties.fragmentTolerance = fragmentTolerance;

        /* check if our directories exist (no typos...) */
        boolean fileDoesNotExist = false;
        for (File spectraDirectoryOrFile : Properties.spectraDirectoryOrFileList) {
            if (!spectraDirectoryOrFile.exists()) {
                U.p("ERROR!  This spectrum file does not exist: " + spectraDirectoryOrFile.getAbsolutePath());
                fileDoesNotExist = true;
            }
        }
        for (File sequenceDirectoryOrFile : Properties.sequenceDirectoryOrFileList) {
            if (!sequenceDirectoryOrFile.exists()) {
                U.p("ERROR!  This sequence file does not exist: " + sequenceDirectoryOrFile.getAbsolutePath());
                fileDoesNotExist = true;
            }
        }
        if (fileDoesNotExist) throw new Error("Search file not found");



        /* we shall use this variable to track how many searches we performed.
         * this will be used in the labeling of our report directory
         */
        int reportIndex = 0;



        /* if there re multiple jobs, the latter varimod settings will persist.
         * we reset those to normal, modification-less parameters
         */
        Properties.scoringMethodName = "Peppy.Match_IMP";
        Properties.matchConstructor = new MatchConstructor(Properties.scoringMethodName);
        Properties.searchModifications = false;

        for (File spectraDirectoryOrFile : Properties.spectraDirectoryOrFileList) {
            Properties.spectraDirectoryOrFile = spectraDirectoryOrFile;
            ArrayList<Spectrum> spectra = SpectrumLoader.loadSpectra();
            int originalSpectraSize = spectra.size();



            /* we are going to combine all returned match sets here */
            ArrayList<Match> allMatchesForSpectralSet = new ArrayList<Match>();

            /* iterate through all of our peptide sources */
            for (int sequenceIndex = 0; sequenceIndex < Properties.sequenceDirectoryOrFileList.size(); sequenceIndex++) {

                /* if we are out of spectra (not likely), get out of this loop */
                if (spectra.size() == 0) break;

                /* set up our sequence data.
                 * Setting this property will also affect the FDR calculation */
                Properties.sequenceDirectoryOrFile = Properties.sequenceDirectoryOrFileList.get(sequenceIndex);
                Properties.isSequenceFileNucleotide = Properties.isSequenceFileNucleotideList.get(sequenceIndex);
                ArrayList<Sequence> sequences = Sequence.loadSequenceFiles(Properties.sequenceDirectoryOrFile);

                FDR fdr = new FDR(spectra);



                /* Calculate score threshold with FDR.
                 * maximumFDR may be negative, indicating we won't use FDR to calculate score cutoff */
                if (Properties.maximumFDR >= 0) {
                    double potentialMinimumScore = fdr.getScoreThreshold(Properties.maximumFDR);

                    /* if we will not find any matches with confidence, skip this round */
                    //NOTE:  in the event of "continue" it will produce no report.  Look out for this when assembling reports!
                    if (potentialMinimumScore < 0) continue;
                    Properties.minimumScore = potentialMinimumScore;
                }



                /* create spectra-based containers for matches */
                ArrayList<MatchesSpectrum> matchesSpectra;
                ArrayList<Match> matches = null;


                matchesSpectra = fdr.getSpectraMatches();
                matches = Peppy.getMatchesFromSpectraMatches(matchesSpectra);


                /* we found no matches */
                if (matches == null) continue;



                /* label the identified peptides according to our sequence */
                for (Match match : matches) {
                    match.getPeptide().setTrackingIdentifier(sequenceIndex);
                }

                /* add these matches to our large pile of all matches for this spectral set */
                allMatchesForSpectralSet.addAll(matches);

                /* group spectra that have been identified */
                Hashtable<Integer, Integer> spectrumIDs = new Hashtable<Integer, Integer>(spectra.size());
                for (Match match : matches) {
                    spectrumIDs.put(match.getSpectrum().getId(), match.getSpectrum().getId());
                }
                U.p(fragmentTolerance + "\t" + spectrumIDs.size());

            }
        }

    }

    public static void runPeppy(String[] args) {

        runFunnelPeppy(args);
    }


    public static void runJobs(String[] args) {
        /* see if we have jobs in the jobs folder */
        File jobsDir = new File("jobs");
        File[] potentialJobsFiles = jobsDir.listFiles();
        ArrayList<File> jobFiles = new ArrayList<File>();
        if (potentialJobsFiles != null) {
            for (int i = 0; i < potentialJobsFiles.length; i++) {
                if (potentialJobsFiles[i].getName().toLowerCase().endsWith(".txt")) {
                    jobFiles.add(potentialJobsFiles[i]);
                }
            }
        }

        /* run the jobs */
        if (jobFiles.size() == 0) {
            U.p("no jobs in jobs folder.  running according to main properties file");
            runPeppy(null);
        } else {
            for (int i = 0; i < jobFiles.size(); i++) {
                propertyFileString = jobFiles.get(i).getName();

                /*
                 * this initializes with the default "properties.txt"
                 * provides a base so common properties don't need to be
                 * redefined
                 */
                Peppy.init(args);

                /*
                 * All properties defined in this file override those defined in "properties.txt"
                 */
                Peppy.init(jobFiles.get(i).getAbsolutePath());

                /* try running Peppy */
                try {
                    runPeppy(null);
                } catch (Exception e) {
                    U.p("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                    e.printStackTrace();
                    U.p("\r\r");
                    U.p("An error occurred. Continuing with the next job...\r\r");
                }
            }
        }
    }

}
