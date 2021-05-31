package Peppy;

import Graphs.Point;
import Math.MassError;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;

public class OptimalTolerances {

    ArrayList<Double> precursorErrors = new ArrayList<Double>();
    ArrayList<Double> fragmentErrors = new ArrayList<Double>();

    ArrayList<Point> precursorScatter = new ArrayList<Point>();

    ArrayList<Double> absolutePrecursorErrors = new ArrayList<Double>();
    ArrayList<Double> absoluteFragmentErrors = new ArrayList<Double>();
    int minimumNumberOfPeaks = Integer.MAX_VALUE;
    double minimumMass = Double.MAX_VALUE;
    double maximumMass = Double.NEGATIVE_INFINITY;
    double optimalPrecursorError;
    double meanPrecursorError;
    double optimalFragmentError;
    double meanFragmentError;

    /* no matches may be found at a desired FDR; this reports that */
    boolean isValid = true;

    private final double portionOfPrecursorErrorsToPass = .99;
    private final double portionOfFragmentErrorsToPass = .99;
    private final double desiredFDR = .01;


    public OptimalTolerances(FDR fdr) {

        double scoreThreshold = fdr.getScoreThreshold(desiredFDR);

        /* if no good confidence at given false discovery rate */
        if (scoreThreshold < 0) {
            isValid = false;
            return;
        }

        ArrayList<MatchesSpectrum> spectraMatches = fdr.getSpectraMatches();

        for (MatchesSpectrum matchesSpectrum : spectraMatches) {

            /* test to see if we should count this datum */
            if (matchesSpectrum.getScore() < scoreThreshold) break;
            if (matchesSpectrum.getMatches().size() == 0) break;
            Match match = matchesSpectrum.getMatches().get(0);
            if (match.getPeptide().isDecoy()) continue;

            /* add to our list of precursor errors */
            double precursorError = MassError.getPPMDifference(match.getPeptide().getMass(), match.getSpectrum().getMass());
            precursorErrors.add(precursorError);
            absolutePrecursorErrors.add(Math.abs(precursorError));

            /* adding to the precursor scatter plot */
            precursorScatter.add(new Point(match.getPeptide().getMass(), precursorError));

            /* harvesting the fragment error sets */
            ArrayList<Double> matchFragmentErrors = match.getFragmentErrors();
            fragmentErrors.addAll(matchFragmentErrors);
            for (double fragmentError : matchFragmentErrors) {
                absoluteFragmentErrors.add(Math.abs(fragmentError));
            }


            if (match.getSpectrum().getPeaks().size() < minimumNumberOfPeaks)
                minimumNumberOfPeaks = match.getSpectrum().getPeaks().size();
            if (match.getPeptide().getMass() < minimumMass) minimumMass = match.getPeptide().getMass();
            if (match.getPeptide().getMass() > maximumMass) maximumMass = match.getPeptide().getMass();
        }

        int cutoffIndex;

        Collections.sort(precursorErrors);
        Collections.sort(fragmentErrors);

        U.p("Min precursor: " + precursorErrors.get(0));
        U.p("max precursor: " + precursorErrors.get(precursorErrors.size() - 1));

        U.p("Min fragment: " + fragmentErrors.get(0));
        U.p("max fragment: " + fragmentErrors.get(fragmentErrors.size() - 1));

        /* finding optimal precursor error */
        Collections.sort(absolutePrecursorErrors);
        cutoffIndex = (int) Math.round((double) absolutePrecursorErrors.size() * portionOfPrecursorErrorsToPass);
        if (cutoffIndex >= absolutePrecursorErrors.size()) cutoffIndex = absolutePrecursorErrors.size() - 1;
        optimalPrecursorError = absolutePrecursorErrors.get(cutoffIndex);



        /* finding optimal fragment error */
        Collections.sort(absoluteFragmentErrors);
        cutoffIndex = (int) Math.round((double) absoluteFragmentErrors.size() * portionOfFragmentErrorsToPass);
        if (cutoffIndex >= absoluteFragmentErrors.size()) cutoffIndex = absoluteFragmentErrors.size() - 1;
        optimalFragmentError = absoluteFragmentErrors.get(cutoffIndex);

        /* finding average fragment error */
        meanFragmentError = 0;
        for (Double fragmentError : fragmentErrors) {
            meanFragmentError += fragmentError;
        }
        meanFragmentError /= fragmentErrors.size();


        /* finding average precursor error */
        meanPrecursorError = 0;
        for (Double precursorError : precursorErrors) {
            meanPrecursorError += precursorError;
        }
        meanPrecursorError /= precursorErrors.size();


    }


    public void createReport(File folder) {

        int numberOfBars = 101;

        double max = Math.max(optimalFragmentError, optimalPrecursorError);
        numberOfBars = 2 * (int) max;



        /* creating the axis */
        double[] precursorAxis = new double[numberOfBars];
        int axisMinimum = -numberOfBars / 2;
        for (int i = 0; i < numberOfBars; i++) {
            precursorAxis[i] = axisMinimum + i;
        }

        /* creating the histograms */
        double[] precursorHistogram = getHistogram(precursorErrors, numberOfBars);
        double[] fragmentHistogram = getHistogram(fragmentErrors, numberOfBars);


        try {
            /* creating the report */
            PrintWriter errorReport = new PrintWriter(new BufferedWriter(new FileWriter(new File(folder, "errorReport.html"))));


            errorReport.println("<html>");
            errorReport.println("  <head>");
            errorReport.println("    <script type=\"text/javascript\" src=\"https://www.google.com/jsapi\"></script>");
            errorReport.println("    <script type=\"text/javascript\">");
            errorReport.println("      google.load(\"visualization\", \"1\", {packages:[\"corechart\"]});");
            errorReport.println("      google.setOnLoadCallback(drawChart);");
            errorReport.println("      function drawChart() {");
            errorReport.println("        var data = google.visualization.arrayToDataTable([");
            errorReport.println("          ['PPM', 'Precursor', 'Fragment'],");

            for (int i = 0; i < numberOfBars; i++) {
                String printString = "[" + precursorAxis[i] + ", " + precursorHistogram[i] + ", " + fragmentHistogram[i] + "]";
                if (i != numberOfBars - 1) {
                    printString += ",";
                }
                errorReport.println(printString);
            }

            errorReport.println("        ]);");
            errorReport.println("");
            errorReport.println("        var options = {");
            errorReport.println("          title: 'Mass Error Distributions',");
            errorReport.println("          hAxis: {title: 'PPM',  titleTextStyle: {color: 'red'}},");
            errorReport.println("          hAxis: {baseline: 0}");
            errorReport.println("        };");
            errorReport.println("");
            errorReport.println("        var chart = new google.visualization.AreaChart(document.getElementById('chart_div'));");
            errorReport.println("        chart.draw(data, options);");
            errorReport.println("      }");
            errorReport.println("    </script>");
            errorReport.println("  </head>");
            errorReport.println("  <body>");
            errorReport.println("    <div id=\"chart_div\" style=\"width: 900px; height: 500px;\"></div>");
            errorReport.println("  </body>");
            errorReport.println("</html>");

            errorReport.flush();
            errorReport.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }


        try {
            /* creating the report */
            PrintWriter errorReport = new PrintWriter(new BufferedWriter(new FileWriter(new File(folder, "precursorScatter.html"))));


            errorReport.println("<html>");
            errorReport.println("  <head>");
            errorReport.println("    <script type=\"text/javascript\" src=\"https://www.google.com/jsapi\"></script>");
            errorReport.println("    <script type=\"text/javascript\">");
            errorReport.println("      google.load(\"visualization\", \"1\", {packages:[\"corechart\"]});");
            errorReport.println("      google.setOnLoadCallback(drawChart);");
            errorReport.println("      function drawChart() {");
            errorReport.println("        var data = google.visualization.arrayToDataTable([");
            errorReport.println("          ['Precursor in Da', 'Error in PPM'],");

            for (int i = 0; i < precursorScatter.size(); i++) {
                Point point = precursorScatter.get(i);
                String printString = "[" + point.x + ", " + point.y + "]";
                if (i != precursorScatter.size() - 1) {
                    printString += ",";
                }
                errorReport.println(printString);
            }

            errorReport.println("        ]);");
            errorReport.println("");
            errorReport.println("        var options = {");
            errorReport.println("          title: 'Precursor Mass vs. Precursor Error',");
            errorReport.println("          hAxis: {title: 'Precursor in Da', minValue: 0, maxValue: 15},");
            errorReport.println("          vAxis: {title: 'Error in PPM', minValue: 0, maxValue: 15},");
            errorReport.println("          pointSize: 1,");
            errorReport.println("          legend: 'none'");
            errorReport.println("        };");
            errorReport.println("");
            errorReport.println("        var chart = new google.visualization.ScatterChart(document.getElementById('chart_div'));");
            errorReport.println("        chart.draw(data, options);");
            errorReport.println("      }");
            errorReport.println("    </script>");
            errorReport.println("  </head>");
            errorReport.println("  <body>");
            errorReport.println("    <div id=\"chart_div\" style=\"width: 900px; height: 500px;\"></div>");
            errorReport.println("  </body>");
            errorReport.println("</html>");

            errorReport.flush();
            errorReport.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }

    /**
     * Assumes histogram block size of 1.
     *
     * @param errors
     * @param numberOfBars
     * @return
     */
    private double[] getHistogram(ArrayList<Double> errors, int numberOfBars) {
        double delta = 1;
        double minimum = -(double) numberOfBars / 2.0;
        double maximum = minimum + (delta * numberOfBars);

        /* creating the histogram */
        Collections.sort(errors);
        double[] errorHistogram = new double[numberOfBars];
        double blockUpperBound = minimum + delta;
        int blockIndex = 0;
        for (double error : errors) {
            if (error < minimum) continue;
            if (error > maximum) continue;
            while (error > blockUpperBound) {
                blockUpperBound += delta;
                blockIndex++;
                if (blockIndex >= numberOfBars) break;
            }
            if (blockIndex >= numberOfBars) break;
            errorHistogram[blockIndex]++;
        }

        /* normalizing the histogram */
        double maxBarValue = 0;
        for (double barValue : errorHistogram) {
            if (barValue > maxBarValue) maxBarValue = barValue;
        }
        for (int i = 0; i < numberOfBars; i++) {
            errorHistogram[i] /= maxBarValue;
        }

        return errorHistogram;
    }


    public int getMinimumNumberOfPeaks() {
        return minimumNumberOfPeaks;
    }


    public double getMinimumMass() {
        return minimumMass;
    }


    public double getMaximumMass() {
        return maximumMass;
    }


    public double getOptimalPrecursorError() {
        return optimalPrecursorError;
    }


    public double getOptimalFragmentError() {
        return optimalFragmentError;
    }


    public boolean isValid() {
        return isValid;
    }


    public double getMeanFragmentError() {
        return meanFragmentError;
    }


    public double getMeanPrecursorError() {
        return meanPrecursorError;
    }

}
