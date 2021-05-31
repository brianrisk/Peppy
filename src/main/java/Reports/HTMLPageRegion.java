package Reports;

import Peppy.Match;
import Peppy.Region;

import java.io.File;
import java.util.ArrayList;

/**
 * Copyright 2013, Brian Risk
 *
 * @author Brian Risk
 */
public class HTMLPageRegion extends HTMLPageMatches {

    Region region;

    public HTMLPageRegion(File destinationFile, Region region) {
        super(region.getMatches(), region.getMatches().size(), destinationFile);
        this.region = region;
    }

    @Override
    public void makePage() {
        String regionName = region.getSequence().getSequenceFile().getName() + ": " + region.getStartLocation() + " to " + region.getStopLocation();
        printHeader(regionName, "<script src=\"http://geneffects.com/files/peppy/js/processing.js\"></script>");
        printP("Sequence: " + region.getSequence().getSequenceFile().getName());

        /* UCSC link */
        printP("region: " + region.getStartLocation() + " to " + region.getStopLocation());
        printP(UCSC.getButton(region.getStartLocation(), region.getStopLocation(), region.getSequence(), matches));

        /* get access to the matches */
        ArrayList<Match> matches = region.getMatches();

        /* print out the processing section */
        int frameHeight = 15;
        int processngHeight = 6 * frameHeight;
        int processingWidth = 683; /* 683 is the approximate width of the UCSC browser */
        print("<script type=\"application/processing\">");
        print("void setup() {");
        print("noLoop();");
        print("size(" + processingWidth + ", " + processngHeight + ");");

        print("}");
        print("void draw() {");
        print("rect(0,0," + (processingWidth - 1) + " ," + (processngHeight - 1) + ");");

        print("noStroke();");

        /* color to mark reverse frame */
        print("fill(255,255,0, 64);"); //a 25% yellow
        print("rect(1,1," + (processingWidth - 1) + "," + (3 * frameHeight - 1) + ");");

        print("fill(0, 64);"); //a 25% black
        int x, y, width;
        y = 0;
        for (Match match : matches) {
            x = match.getPeptide().getStartIndex() - region.getStartLocation();
            x = scaleInt(x, region.getMaxLength(), processingWidth);
            y = (match.getPeptide().getStartIndex() % 3) * frameHeight;
            if (match.getPeptide().isForward()) y += frameHeight * 3;
            width = match.getPeptide().getStopIndex() - match.getPeptide().getStartIndex();
            width = scaleInt(width, region.getMaxLength(), processingWidth);
            print("rect(" + x + ", " + y + ", " + width + ", " + frameHeight + ");");
        }
        /* draw the lines delineating the different frames */
        print("stroke(64);");
        for (int i = 0; i < 5; i++) {
            int lineLevel = (i + 1) * frameHeight;
            print("line(1," + lineLevel + "," + (processingWidth - 1) + "," + lineLevel + ");");
        }
        print("}");
        print("</script><canvas width=\"" + region.getMaxLength() + "px\" height=\"" + processngHeight + "px\"></canvas>");


        printH2("Best Matches");
        printBestMatchTable("../../spectra/");
        printH2("All Matches");
        printMatchTable("../../spectra/");


        printFooter();

    }

    /**
     * math func to let us scale coordinates when we are drawing a region smaller (or larger) than it is.
     *
     * @param x
     * @param oldWidth
     * @param newWidth
     * @return
     */
    private int scaleInt(int x, int oldWidth, int newWidth) {
        double dx = x;
        dx *= newWidth;
        dx /= oldWidth;
        return (int) Math.round(dx);
    }

}
