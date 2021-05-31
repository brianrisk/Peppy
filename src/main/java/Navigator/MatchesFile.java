package Navigator;

import java.io.File;

/**
 * Copyright 2013, Brian Risk
 *
 * @author Brian Risk
 */
public class MatchesFile {

    File file;


    private static int typeTracker = 0;
    public static final int TYPE_PROTEIN = typeTracker++;
    public static final int TYPE_DNA = typeTracker++;
    public static final int TYPE_RNA = typeTracker++;

}
