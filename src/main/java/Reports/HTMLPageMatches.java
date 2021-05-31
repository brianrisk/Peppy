package Reports;

import Peppy.Match;
import Peppy.Matches;
import Peppy.Properties;

import java.io.File;
import java.text.NumberFormat;
import java.util.ArrayList;

/**
 * Copyright 2013, Brian Risk
 *
 * @author Brian Risk
 */
public class HTMLPageMatches extends HTMLPage {

    ArrayList<Match> matches;
    int maxDisplay;

    public HTMLPageMatches(ArrayList<Match> matches, int maxDisplay, File destinationFile) {
        super(destinationFile);
        this.matches = matches;
        this.maxDisplay = maxDisplay;
    }

    @Override
    public void makePage() {

        printHeader();

        printLink("regions html/index.html", "View regions report");

        /* print table of best matches */
        printH2("Best matches for each peptide");
        printBestMatchTable("spectra/");

        printFooter();
    }

    protected void printBestMatchTable(String spectraPath) {
        printMatchTable(spectraPath, Matches.getBestMatches(matches));
    }

    protected void printMatchTable(String spectraPath) {
        printMatchTable(spectraPath, matches);
    }

    protected void printMatchTable(String spectraPath, ArrayList<Match> matches) {
        /* print our table headers */
        print("<table class=\"sortable\" id=\"box-table-a\" width=\"95%\">");
        printTR();
        printTH("ID");
        printTH("acid");
        printTH("location");
        printTH("score");
        printTH("has mod");
        printTH("count");

        /*print the rows */
        printMatchRows(matches, spectraPath);

        print("</table>");
    }


    protected void printMatchRows(ArrayList<Match> matches, String spectraPath) {
        /* some basic bounds checks */
        if (matches == null) return;
        if (matches.size() == 0) return;
        int maxDisplay = this.maxDisplay;
        if (matches.size() < maxDisplay) maxDisplay = matches.size();

        /* set up how we will round */
        NumberFormat nfDecimal = NumberFormat.getInstance();
        nfDecimal.setMaximumFractionDigits(2);

        for (int i = 0; i < maxDisplay; i++) {
            /* get our match */
            Match match = matches.get(i);

            printTR();

            /* spectrum ID link */
            String spectrumLink = "<a href=\"" + spectraPath + match.getSpectrum().getId() + Properties.reportWebSuffix + "\">" + match.getSpectrum().getId() + "</a>";
            printTD(spectrumLink);

            /* the acid string */
            String peptideMarkup = "";
            peptideMarkup += match.getPeptide().getAcidSequenceString();
            /* google link */
            String link = "http://www.google.com/search?hl=en&q=" + match.getPeptide().getAcidSequenceString();
            printTD("(<a href=\"" + link + "\">" + peptideMarkup + "</a>)");

            /* score */
            printTD(nfDecimal.format(match.getScore()));

            /* modifications */
            printTD("" + match.hasModification());

            /* count */
            printTD("" + match.getSpectrumMatches().getMatches().size());
        }
    }

}
