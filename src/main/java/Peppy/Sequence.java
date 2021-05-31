package Peppy;

import java.io.File;
import java.util.ArrayList;

/**
 * Copyright 2013, Brian Risk
 *
 * @author Brian Risk
 */
public abstract class Sequence {

    protected File sequenceFile;
    protected int id;


    /**
     * The amount of peptides produced by a sequence may be huge.
     * therefore, this method is set up so that it should be called
     * multiple times, each time it is called it returns another portion of the
     * database of peptides.  When there are no more peptides from that sequence
     * then null is returned.
     *
     * @param reverse This is if we want to reverse our database to get a null database
     * @return a sorted ArrayList of amino acid sequence fragments from the given sequence file
     */
    public abstract ArrayList<Peptide> extractMorePeptides(boolean reverse);

    /**
     * Runs extractMorePeptides() until there are no more peptides
     * to extract.  Returns the full list.
     *
     * @param reverse This is if we want to reverse our database to get a null database
     * @return
     */
    public abstract ArrayList<Peptide> extractAllPeptides(boolean reverse);

    public abstract void reset();


    public File getSequenceFile() {
        return sequenceFile;
    }

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    public static ArrayList<Sequence> loadSequenceFiles(File folder) {
        ArrayList<Sequence> sequences = new ArrayList<Sequence>();
        if (folder.isFile()) {
            if (Properties.isSequenceFileNucleotide) {
                sequences.add(new SequenceNucleotide(folder));
            } else {
                sequences.add(new SequenceProtein(folder));
            }
        } else {
            File[] files = folder.listFiles();
            for (int i = 0; i < files.length; i++) {
                if (files[i].isHidden()) continue;
                String fileName = files[i].getName().toLowerCase();
                if (
                        fileName.endsWith(".fasta") ||
                                fileName.endsWith(".fa") ||
                                fileName.endsWith(".fsa") ||
                                fileName.endsWith(".dat") ||
                                fileName.endsWith(".txt")) {
                    if (Properties.isSequenceFileNucleotide) {
                        sequences.add(new SequenceNucleotide(files[i]));
                    } else {
                        sequences.add(new SequenceProtein(files[i]));
                    }
                }
            }
        }
        for (int i = 0; i < sequences.size(); i++) {
            sequences.get(i).setId(i);
        }
        return sequences;
    }

    public static ArrayList<Sequence> loadSequenceFiles(String folderName) {
        return loadSequenceFiles(new File(folderName));
    }


}
