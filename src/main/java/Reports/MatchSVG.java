package Reports;

import Math.MassError;
import Peppy.*;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

public class MatchSVG {

    Match match;
    File outFile;
    PrintWriter pw;
    Spectrum spectrum;
    Peptide peptide;
    double massRange;
    double spectrumMaxIntensity;
    String acidSequence;
    ArrayList<PeakLabeled> labeledPeaks = new ArrayList<PeakLabeled>();

    //Image specs
    int chartWidth = 400;
    int chartHeight = 300;
    int chartX = 45;
    int chartY = 20;
    int bottomPadding = 35;
    int spectrumWidth = chartWidth - chartX;
    int spectrumHeight = chartHeight - chartY - bottomPadding;

    //track if we need to redraw the image;
    boolean redrawRequired = true;

    //font properties
    int fontHeight = 20;


    boolean displayMasses = false;

    //colors
    String yIonColor = "#FF0000";
    String bIonColor = "#0000FF";
    String noIonColor = "#888888";
    String bothIonColor = "#00FF00";
    String bkgndColor = "#FFFFFF";
    String axesColor = "#000000";

    String strokeColor = "#000000";
    String fillColor = "#000000";
    double strokeWeight = 1;

    int fontSize = 9;
    Font font = new Font("Helvetica", Font.PLAIN, fontSize);
    FontMetrics fontMetrics;

    public static void main(String[] args) {
        Match_Blank match = new Match_Blank(new Spectrum("/Users/risk2/Documents/Papers/Peppy Hourglass method paper/COBRA1 ASGVSAAAPGER  peptide/Research_302_2984_261_69_2_columns.0718.0718.2.dta"), new Peptide("ASGVSAAAPGER"), 50);
        MatchSVG msvg = new MatchSVG(match, new File("ASGVSAAAPGER.svg"));
        msvg.saveSVG();
        U.p("done");
    }

    public MatchSVG(Match match, File outFile) {
        this.match = match;
        try {
            pw = new PrintWriter(new FileWriter(outFile));
        } catch (IOException e) {
            e.printStackTrace();
        }
        pw.println("<svg xmlns=\"http://www.w3.org/2000/svg\">");
        this.spectrum = match.getSpectrum();
        this.peptide = match.getPeptide();
        massRange = spectrum.getMass();
        spectrumMaxIntensity = spectrum.getMaxIntensity();
        acidSequence = peptide.getAcidSequenceString();


        int chartWidth = spectrumWidth + 100;
        int chartHeigth = spectrumHeight + 50;
        BufferedImage chartImage = new BufferedImage(chartWidth, chartHeigth, BufferedImage.TYPE_INT_RGB);
        Graphics2D chartGraphics = chartImage.createGraphics();
        fontMetrics = chartGraphics.getFontMetrics(font);

        ArrayList<Peak> peaks = spectrum.getPeaks();
        for (Peak peak : peaks) {
            labeledPeaks.add(new PeakLabeled(peak));
        }

        markMatchingIons();

        /*
         * Bit is a hack to sort by intensity for the Peppy figure
         */
//		for (Peak peak: labeledPeaks) {
//			peak.setCompareByIntensity();
//		}
//		Collections.sort(labeledPeaks);
//		
//		for (int i = 0; i < labeledPeaks.size(); i++) {;
//			PeakLabeled peak = labeledPeaks.get(i);
//			peak.setMass(i * 12);
//		}


    }

    private void line(double x1, double y1, double x2, double y2) {
        pw.println("<line x1=\"" + x1 + "\" y1=\"" + y1 + "\" x2=\"" + x2 + "\" y2=\"" + y2 + "\" stroke-width=\"" + strokeWeight + "\" stroke=\"" + strokeColor + "\" />");
    }

    private void rect(double x, double y, double width, double height) {
        pw.println("<rect x=\"" + x + "\" y=\"" + y + "\" width=\"" + width + "\" height=\"" + height + "\" stroke-width=\"" + strokeWeight + "\" stroke=\"" + strokeColor + "\" fill=\"" + fillColor + "\" />");
    }

    private void text(String text, double x, double y) {
        pw.println("<text x=\"" + x + "\" y=\"" + y + "\" font-size=\"" + fontSize + "\" fill=\"" + fillColor + "\">" + text + "</text>");
    }

    private void stroke(String color) {
        strokeColor = color;
    }

    private void strokeWeight(double x) {
        strokeWeight = x;
    }

    private void fill(String color) {
        fillColor = color;
    }

    private double log(double number) {
        return Math.log(number);
    }

    private double round(double x) {
        return Math.round(x);
    }

    private double pow(double x, double exp) {
        return Math.pow(x, exp);
    }

    private int textWidth(String string) {
        return fontMetrics.stringWidth(string);
    }

    public void saveSVG() {

        //draw the spectrum
        drawSpectrum();

        //draw axes
        stroke(axesColor);
        line(chartX, chartY + spectrumHeight + 1, chartX + spectrumWidth, chartY + spectrumHeight + 1); //x axis
        line(chartX + spectrumWidth, chartY + spectrumHeight - 10, chartX + spectrumWidth, chartY + spectrumHeight + 10); //x axis cap
        line(chartX, chartY, chartX, chartY + spectrumHeight + 2); //y axis
        line(chartX - 10, chartY, chartX + 10, chartY); // y cap

        //find x tick marks
        double xLabelIncrement = log(massRange);
        xLabelIncrement /= log(10);
        xLabelIncrement = pow(10, round(xLabelIncrement) - 1);
        int xPixelIncrement = (int) ((xLabelIncrement / massRange) * spectrumWidth);
        //to avoid crowding
        if (xPixelIncrement < 25) {
            xLabelIncrement = log(massRange);
            xLabelIncrement /= log(10);
            xLabelIncrement = pow(10, round(xLabelIncrement));
            xPixelIncrement = (int) ((xLabelIncrement / massRange) * spectrumWidth);
        }

        //draw x ticks
        int x1;
        int numberOfTicks = (spectrumWidth / xPixelIncrement);
        fill("#440088");
        for (int i = 0; i <= numberOfTicks; i++) {
            x1 = chartX + i * xPixelIncrement;
            line(x1, chartY + spectrumHeight + 1, x1, chartY + spectrumHeight + 12);
            String tickLabel = "" + ((int) (xLabelIncrement * i));
            text(tickLabel, x1 - 5, chartY + spectrumHeight + fontHeight + 0);
        }

        //find y tick marks
        double yLabelIncrement = log(spectrumMaxIntensity);
        yLabelIncrement /= log(10);
        yLabelIncrement = pow(10, round(yLabelIncrement) - 1);
        int yPixelIncrement = (int) ((yLabelIncrement / spectrumMaxIntensity) * spectrumHeight);
        //to avoid crowding
        if (yPixelIncrement < 10) {
            yLabelIncrement = log(spectrumMaxIntensity);
            yLabelIncrement /= log(10);
            yLabelIncrement = pow(10, round(yLabelIncrement));
            yPixelIncrement = (int) ((yLabelIncrement / spectrumMaxIntensity) * spectrumHeight);
        }

        //draw y ticks
        int y1;
        numberOfTicks = (spectrumHeight / yPixelIncrement);

        for (int i = 1; i <= numberOfTicks; i++) {
            double labelValue = yLabelIncrement * i;
            y1 = chartY + spectrumHeight - i * yPixelIncrement;
            line(chartX - 10, y1, chartX, y1);
            String tickLabel = "" + ((int) (labelValue));
            int labelWidth = textWidth(tickLabel);
            text(tickLabel, chartX - 12 - labelWidth, y1);
        }

        //draw labels
        fill("#FF0000");
        text(acidSequence, chartX + 10, chartY + fontHeight);

        pw.println("</svg>");
        pw.flush();
        pw.close();

    }


    /*
     * this should be how short it is in Java!
     */
    void drawSpectrum() {
        double xLoc, yLoc;

        //getting maximum spectrum value and intensity
        int maxValue = (int) massRange;
        int massStart = 0;
        int massStop = maxValue;

        /* non-matches layer */
        strokeWeight(1);
        for (PeakLabeled peak : labeledPeaks) {
            if (peak.getMass() >= massStart && peak.getMass() <= massStop) {
                if (peak.getColor() == noIonColor) {
                    stroke(peak.getColor());
                    xLoc = ((peak.getMass() - massStart) * spectrumWidth / massRange);
                    double dropAmount = 0;
                    dropAmount = ((peak.getIntensity() / spectrumMaxIntensity) * spectrumHeight);
                    if (dropAmount < 2) dropAmount = 2;
                    yLoc = (spectrumHeight - dropAmount);


                    line(xLoc + chartX, yLoc + chartY, xLoc + chartX, spectrumHeight + chartY);
                }
            }
        }

        /* matches layer */
        strokeWeight(1.5);
        for (PeakLabeled peak : labeledPeaks) {
            if (peak.getMass() >= massStart && peak.getMass() <= massStop) {
                if (peak.getColor() != noIonColor) {
                    stroke(peak.getColor());
                    xLoc = ((peak.getMass() - massStart) * spectrumWidth / massRange);
                    double dropAmount = 0;
                    dropAmount = ((peak.getIntensity() / spectrumMaxIntensity) * spectrumHeight);
                    if (dropAmount < 2) dropAmount = 2;
                    yLoc = (spectrumHeight - dropAmount);


                    strokeWeight(1.5);
                    if (yLoc + chartY < fontHeight) yLoc += fontHeight;
                    String label = "";
                    if (peak.getbIonNumber() != 0) label += "b" + peak.getbIonNumber();
                    if (peak.getbIonNumber() != 0 && peak.getyIonNumber() != 0) label += ", ";
                    if (peak.getyIonNumber() != 0) label += "y" + peak.getyIonNumber();
                    if (displayMasses) label += ", " + peak.getMass();
                    fill(peak.getColor());
                    text(label, xLoc + 4 + chartX, yLoc + chartY);

                    line(xLoc + chartX, yLoc + chartY, xLoc + chartX, spectrumHeight + chartY);
                }
            }
        }
    }


    /*
     * Used for the figure in the peppy score paper
     * explaining the firs, random alignment p-value.
     */
    void drawSpectrumWithGreyBoxes() {
        double xLoc, yLoc;

        //getting maximum spectrum value and intensity
        int maxValue = (int) massRange;
        int massStart = 0;
        int massStop = maxValue;

        /* non-mathches layer */
        strokeWeight(0);
        for (PeakLabeled peak : labeledPeaks) {
            if (peak.getMass() >= massStart && peak.getMass() <= massStop) {
                xLoc = ((peak.getMass() - massStart) * spectrumWidth / massRange);
                double dropAmount = 0;
                dropAmount = ((peak.getIntensity() / spectrumMaxIntensity) * spectrumHeight);
                if (dropAmount < 2) dropAmount = 2;
                yLoc = (spectrumHeight - dropAmount);

                fill("#CCCCCC");
                rect(xLoc + chartX - 2, chartY, 4, spectrumHeight);
            }
        }
        strokeWeight(1);
        for (PeakLabeled peak : labeledPeaks) {
            if (peak.getMass() >= massStart && peak.getMass() <= massStop) {
                xLoc = ((peak.getMass() - massStart) * spectrumWidth / massRange);
                double dropAmount = 0;
                dropAmount = ((peak.getIntensity() / spectrumMaxIntensity) * spectrumHeight);
                if (dropAmount < 2) dropAmount = 2;
                yLoc = (spectrumHeight - dropAmount);

                stroke("#000000");
                line(xLoc + chartX, yLoc + chartY, xLoc + chartX, spectrumHeight + chartY);
            }
        }

    }


    /* mirror these changes in Peppy */
    private void markMatchingIons() {
        int i;
        double theoreticalPeakMass, peakMass;
        int peakIndex, seqIndex;

        //initialize peakColors
        for (PeakLabeled peak : labeledPeaks) {
            peak.setColor(noIonColor);
        }
        if (acidSequence.length() == 0) return;


        //we want -1 because most of these spectra will have a match with
        //the last theoretical peak
        int peptideLengthMinusOne = acidSequence.length() - 1;
        if (acidSequence.charAt(peptideLengthMinusOne) == '.') peptideLengthMinusOne--;

        //find the ranges around our theoretical peptides where we
        //count spectrum peaks
        double[] theoreticalPeaksLeft = new double[peptideLengthMinusOne];
        double[] theoreticalPeaksRight = new double[peptideLengthMinusOne];

        //y ion
        //computing the left and right boundaries for the ranges where our peaks should land
        //var theoreticalPeakMass = getPeptideMass() + rightIonDifference;
        theoreticalPeakMass = peptide.getMass() + Properties.rightIonDifference;


        for (i = 0; i < peptideLengthMinusOne; i++) {
            theoreticalPeakMass -= peptide.getResidueMass(i);
            theoreticalPeaksLeft[i] = theoreticalPeakMass - MassError.getDaltonError(Properties.fragmentTolerance, theoreticalPeakMass);
            theoreticalPeaksRight[i] = theoreticalPeakMass + MassError.getDaltonError(Properties.fragmentTolerance, theoreticalPeakMass);
        }

        peakIndex = labeledPeaks.size() - 1;
        seqIndex = 0;
        while (peakIndex >= 0) {
            PeakLabeled peak = labeledPeaks.get(peakIndex);
            peakMass = peak.getMass();

            while (peakMass < theoreticalPeaksLeft[seqIndex]) {
                seqIndex++;
                if (seqIndex == peptideLengthMinusOne) break;
            }

            if (seqIndex == peptideLengthMinusOne) break;
            if (peakMass >= theoreticalPeaksLeft[seqIndex] && peakMass <= theoreticalPeaksRight[seqIndex]) {
                peak.setColor(yIonColor);
                peak.setyIonNumber(peptideLengthMinusOne - seqIndex);
            }

            peakIndex--;
        }

        //b ion
        theoreticalPeakMass = Properties.leftIonDifference;
        for (i = 0; i < peptideLengthMinusOne; i++) {
            theoreticalPeakMass += peptide.getResidueMass(i);
            theoreticalPeaksLeft[i] = theoreticalPeakMass - MassError.getDaltonError(Properties.fragmentTolerance, theoreticalPeakMass);
            theoreticalPeaksRight[i] = theoreticalPeakMass + MassError.getDaltonError(Properties.fragmentTolerance, theoreticalPeakMass);
        }

        peakIndex = 0;
        seqIndex = 0;
        while (peakIndex < labeledPeaks.size()) {
            PeakLabeled peak = labeledPeaks.get(peakIndex);
            peakMass = peak.getMass();
            while (peakMass > theoreticalPeaksRight[seqIndex]) {
                seqIndex++;
                if (seqIndex == peptideLengthMinusOne) break;
            }
            if (seqIndex == peptideLengthMinusOne) break;
            if (peakMass >= theoreticalPeaksLeft[seqIndex] && peakMass <= theoreticalPeaksRight[seqIndex]) {
                if (peak.getColor().equals(yIonColor)) {
                    peak.setColor(bothIonColor);
                } else {
                    peak.setColor(bIonColor);
                }
                peak.setbIonNumber(seqIndex + 1);
            }

            peakIndex++;
        }
    }


}
