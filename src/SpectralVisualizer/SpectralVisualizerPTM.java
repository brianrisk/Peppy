package SpectralVisualizer;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.Polygon;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import javax.imageio.ImageIO;

import Peppy.AminoAcids;
import Peppy.Peak;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.Spectrum;

public class SpectralVisualizerPTM  {
	

	private static boolean cubedScale = false;
	static Color yIonColorPRE = Color.RED;
	static Color yIonColorPOST = Color.MAGENTA;
	static Color bIonColorPRE = Color.BLUE;
	static Color bIonColorPOST = Color.CYAN;
	static Color noIonColor = Color.gray;
	static Color bothIonColor = Color.green;
	static Color bkgndColor = Color.white;
	static Color axesColor = Color.black;
	

	
	public static void drawDeluxSpectrum(Spectrum spectrum, Peptide peptide, File imageFile, double offset, int modifiedIndex)  throws IOException {
		//mark matching ions
		markMatchingIons(spectrum, peptide, offset, modifiedIndex);
		
		//getting sectrum image
		int spectrumWidth = 1000;
		int spectrumHeight = 300;
		BufferedImage spectrumImage = new BufferedImage(spectrumWidth, spectrumWidth, BufferedImage.TYPE_INT_RGB);
		drawSpectrum(spectrum, spectrumWidth, spectrumHeight, null, true, "", "", spectrumImage);
		
		//getting full chart graphics
		int chartX = 75;
		int chartY = 10;
		int chartWidth = spectrumWidth + 100;
		int chartHeigth = spectrumHeight + 50;
		BufferedImage chartImage = new BufferedImage(chartWidth, chartHeigth, BufferedImage.TYPE_INT_RGB);
		Graphics2D chartGraphics = chartImage.createGraphics();
		chartGraphics.setBackground(bkgndColor);
		FontMetrics fontMetrics = chartGraphics.getFontMetrics();
		int fontHeight = fontMetrics.getHeight();
		
		//putting on the spectrum
		chartGraphics.setColor(bkgndColor);
		chartGraphics.fillRect(0,0,chartWidth,chartHeigth);
		chartGraphics.drawImage(spectrumImage, chartX, chartY, null);
		
		//fill in a stupid black rectangle
		chartGraphics.fillRect(chartX,chartY + spectrumHeight + 3, spectrumWidth,chartHeigth - chartY - spectrumHeight);
		
		//draw axes
		chartGraphics.setColor(axesColor);
		chartGraphics.drawLine(chartX, chartY + spectrumHeight + 2, chartX + spectrumWidth, chartY + spectrumHeight + 2); //x axis
		chartGraphics.drawLine(chartX + spectrumWidth, chartY + spectrumHeight -10, chartX + spectrumWidth, chartY + spectrumHeight +10); //x axis cap
		chartGraphics.drawLine(chartX, chartY, chartX, chartY + spectrumHeight); //y axis
		chartGraphics.drawLine(chartX - 10, chartY, chartX + 10, chartY); // y cap
		
		//find x tick marks
		double xLabelIncrement = Math.log10(spectrum.getMass());
		xLabelIncrement = Math.pow(10, Math.round(xLabelIncrement) - 1);
		int xPixelIncrement = (int) ((xLabelIncrement / spectrum.getMass()) * spectrumWidth);
		
		//draw x ticks
		int x1;
		int numberOfTicks = (spectrumWidth / xPixelIncrement);
		for (int i = 0; i <= numberOfTicks; i++) {
			x1 = chartX + i * xPixelIncrement;
			chartGraphics.drawLine(x1, chartY + spectrumHeight + 2, x1, chartY + spectrumHeight + 12);
			String tickLabel = "" + (int) xLabelIncrement * i;
			chartGraphics.drawString(tickLabel, x1 - 5, chartY + spectrumHeight + fontHeight + 9);
		}
		
		//find y tick marks
		double maxIntensity = spectrum.getMaxIntensity();
		double yLabelIncrement = Math.log10(maxIntensity);
		yLabelIncrement = Math.pow(10, Math.round(yLabelIncrement) - 1);
		int yPixelIncrement = (int) ((yLabelIncrement / maxIntensity) * spectrumHeight);
		
		//draw y ticks
		int y1;
		numberOfTicks = (spectrumHeight / yPixelIncrement);
		for (int i = 1; i <= numberOfTicks; i++) {
			if (cubedScale) {
				y1 = chartY + spectrumHeight - (int) ((Math.pow(i * yLabelIncrement, 0.333) / Math.pow(maxIntensity, 0.33)) * spectrumHeight);
			} else {
				y1 = chartY + spectrumHeight - i * yPixelIncrement;
			}	
			
			chartGraphics.drawLine(chartX -10, y1, chartX, y1);
			String tickLabel = "" + (int) yLabelIncrement * i;
			int labelWidth = fontMetrics.stringWidth(tickLabel);
			chartGraphics.drawString(tickLabel,chartX - 12 - labelWidth, y1);
		}
		
		//draw labels
		chartGraphics.setColor(Color.RED);
		chartGraphics.drawString(peptide.getAcidSequenceString(), chartX + 10, chartY + fontHeight);
		chartGraphics.drawString(spectrum.getFile().getName(), chartX + 10, (int) (chartY + fontHeight * 2.5));
		
		//save the file
		ImageIO.write(chartImage,"JPG",imageFile);
		
	}
	

	
	public static void drawSpectrum(Spectrum spectrum, int width, int height, File dest, boolean drawLabels, String message1, String message2, BufferedImage bdest) throws IOException {
		
		int xLoc, yLoc;
		Graphics2D g = bdest.createGraphics();
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_SPEED));
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
		g.setColor(bkgndColor);
		g.fillRect(0,0,width,height);
		
		Font oldFont = g.getFont();
		g.setFont(new Font(oldFont.getName(), oldFont.getSize(), (int) (oldFont.getSize() * 1.5)));
		FontMetrics fontMetrics = g.getFontMetrics();
		int fontHeight = fontMetrics.getHeight();
		
		//getting maximum spectrum value and intensity
		double maxValue = spectrum.getMass();
//		double maxValue = 1500;
		double maxIntensity = spectrum.getMaxIntensity();
		if (cubedScale) maxIntensity = Math.pow(maxIntensity, 0.333);
		double peakIntensity;
		
		//first draw the non-important lines
		ArrayList<Peak> peaks = spectrum.getPeaks();
		g.setColor(noIonColor);
		for (Peak peak: peaks) {
			if (peak.getColor().equals(noIonColor)) {
				peakIntensity = peak.getIntensity();
				if (cubedScale) peakIntensity = Math.pow(peakIntensity, 0.333);
				xLoc = (int) (peak.getMass()  * width / maxValue);
				yLoc = (int) (height - ((peakIntensity / maxIntensity) * height));
				g.setStroke(new BasicStroke(1.5f));
				g.drawLine(xLoc, yLoc, xLoc, height);
				
				//if peak is hilighted, draw a little triangle above it
				if (peak.isHilighted()) {
					if (yLoc < 10) yLoc += 10;
					g.setColor(Color.yellow);
					Polygon polygon = new Polygon();
					polygon.addPoint(xLoc - 10, yLoc - 10);
					polygon.addPoint(xLoc - 10, yLoc + 10);
					polygon.addPoint(xLoc, yLoc);
					g.fillPolygon(polygon);
				}
			}
		}
		
		//draw the colored peaks
		for (Peak peak: peaks) {
			if (!peak.getColor().equals(noIonColor)) {
				g.setColor(peak.getColor());
				peakIntensity = peak.getIntensity();
				if (cubedScale) peakIntensity = Math.pow(peakIntensity, 0.333);
				xLoc = (int) (peak.getMass()  * width / maxValue);
				yLoc = (int) (height - ((peakIntensity / maxIntensity) * height));
				g.setStroke(new BasicStroke(2.0f));
				g.drawLine(xLoc, yLoc, xLoc, height);
				if (drawLabels && !peak.getColor().equals(noIonColor)) {
					if (yLoc < fontHeight) yLoc += fontHeight;
					StringBuffer label = new StringBuffer();
					if (peak.getbIonNumber() != -1) label.append("b" + peak.getbIonNumber());
					if (peak.getbIonNumber() != -1 && peak.getyIonNumber() != -1) label.append(", ");
					if (peak.getyIonNumber() != -1) label.append("y" + peak.getyIonNumber());
					g.drawString("" + label, xLoc + 4, yLoc);
				}
				//if peak is hilighted, draw a little triangle above it
				if (peak.isHilighted()) {
					if (yLoc < 10) yLoc += 10;
					g.setColor(Color.yellow);
					Polygon polygon = new Polygon();
					polygon.addPoint(xLoc - 10, yLoc - 10);
					polygon.addPoint(xLoc - 10, yLoc + 10);
					polygon.addPoint(xLoc, yLoc);
					g.fillPolygon(polygon);
				}
			}
		}
		
		//draw messages
		g.setColor(Color.lightGray);
		g.drawString(message1, 10, 20);
		g.drawString(message2, 10, 40);
		
		//save the image
		if (dest != null) {
			ImageIO.write(bdest,"JPG",dest);
		}
	}
	

	
	public static void markMatchingIons(Spectrum spectrum, Peptide peptide, double offset, int modifiedIndex) {
		byte [] acidSequence = peptide.getAcidSequence();
		
		//initialize all peaks to be gray
		ArrayList<Peak> peaks = spectrum.getPeaks();
		for (Peak peak: peaks) {
			peak.setColor(noIonColor);
		}

		int i;
		double theoreticalPeakMass, peakMass;
		int peakIndex, seqIndex;
		
		//we want -1 because most of these spectra will have a match with 
		//the last theoretical peak
		int peptideLengthMinusOne = acidSequence.length - 1;
		
		double [] bIonMatchesWithHighestIntensity = new double[acidSequence.length];
		for (i = 0; i < acidSequence.length; i++) bIonMatchesWithHighestIntensity[i] = 0.0;
		double [] yIonMatchesWithHighestIntensity = new double[acidSequence.length];
		for (i = 0; i < acidSequence.length; i++) yIonMatchesWithHighestIntensity[i] = 0.0;

		
		//find the ranges around our theoretical peptides where we
		//count spectrum peaks
		double [] theoreticalPeaksLeft = new double[peptideLengthMinusOne];
		double [] theoreticalPeaksRight = new double[peptideLengthMinusOne];
		
		
		/* y-ion  */
		//computing the left and right boundaries for the ranges where our peaks should land
		theoreticalPeakMass = offset + peptide.getMass() + Properties.rightIonDifference;
		for (i = 0; i < peptideLengthMinusOne; i++) {
			theoreticalPeakMass -= AminoAcids.getWeightMono(acidSequence[i]);
			if (i == modifiedIndex) theoreticalPeakMass -= offset;
			theoreticalPeaksLeft[i] = theoreticalPeakMass - Properties.peakDifferenceThreshold;
			theoreticalPeaksRight[i] = theoreticalPeakMass + Properties.peakDifferenceThreshold;	
		}
		
		peakIndex = spectrum.getPeakCount() - 1;
		seqIndex = 0;
		while (peakIndex >= 0) {
			peakMass = spectrum.getPeak(peakIndex).getMass();
			while (peakMass < theoreticalPeaksLeft[seqIndex]) {
				seqIndex++;
				if (seqIndex == peptideLengthMinusOne) break;
			}
			if (seqIndex == peptideLengthMinusOne) break;
			if (peakMass >= theoreticalPeaksLeft[seqIndex] && peakMass <= theoreticalPeaksRight[seqIndex]) {
				if (yIonMatchesWithHighestIntensity[seqIndex] < spectrum.getPeak(peakIndex).getIntensity()) {
					spectrum.getPeak(peakIndex).setColor(yIonColorPRE);
					int labelNumber = peptideLengthMinusOne - seqIndex;
//					if (acidSequence.endsWith(".")) labelNumber--;
					spectrum.getPeak(peakIndex).setyIonNumber(labelNumber);
				}
			}
			
			peakIndex--;
		}
			
		
		/* b-ion  */
		theoreticalPeakMass = Properties.leftIonDifference;
		for (i = 0; i < peptideLengthMinusOne; i++) {
			theoreticalPeakMass += AminoAcids.getWeightMono(acidSequence[i]);
			if (i == modifiedIndex) theoreticalPeakMass += offset;
			theoreticalPeaksLeft[i] = theoreticalPeakMass - Properties.peakDifferenceThreshold;
			theoreticalPeaksRight[i] = theoreticalPeakMass + Properties.peakDifferenceThreshold;
		}
		
		peakIndex = 0;
		seqIndex = 0;
		while (peakIndex < spectrum.getPeakCount()) {
			peakMass = spectrum.getPeak(peakIndex).getMass();
			while (peakMass > theoreticalPeaksRight[seqIndex]) {
				seqIndex++;
				if (seqIndex == peptideLengthMinusOne) break;
			}
			if (seqIndex == peptideLengthMinusOne) break;
			if (peakMass >= theoreticalPeaksLeft[seqIndex] && peakMass <= theoreticalPeaksRight[seqIndex]) {
				if (bIonMatchesWithHighestIntensity[seqIndex] < spectrum.getPeak(peakIndex).getIntensity()) {
					if (spectrum.getPeak(peakIndex).getColor().equals(yIonColorPRE)) {
						spectrum.getPeak(peakIndex).setColor(bothIonColor);
					} else {
						spectrum.getPeak(peakIndex).setColor(bIonColorPRE);
					}
					spectrum.getPeak(peakIndex).setbIonNumber(seqIndex + 1);
				}
			}
			
			peakIndex++;
		}
		
	}
	

}
