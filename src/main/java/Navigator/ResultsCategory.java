package Navigator;

import java.io.File;
import java.util.ArrayList;

/**
 * Copyright 2013, Brian Risk
 *
 * @author Brian Risk
 */
public class ResultsCategory {


    int ID;
    String name;
    String description;
    ArrayList<File> files = new ArrayList<File>();
    int databaseType;

    public static final int PROTEIN = 0;
    public static final int RNA = 1;
    public static final int DNA = 2;

    public ResultsCategory(String name, int databaseType) {
        this.name = name;
        this.databaseType = databaseType;
    }

    public String getName() {
        return name;
    }

    public boolean equals(ResultsCategory other) {
        return (name.equals(other.getName()));
    }

    public void addFile(File file) {
        files.add(file);
    }

    public ArrayList<File> getFiles() {
        return files;
    }

    public int getDatabaseType() {
        return databaseType;
    }


}
