package SpectralVisualizer;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.Polygon;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import javax.imageio.ImageIO;

import Peppy.AminoAcids;
import Peppy.Definitions;
import Peppy.Peak;
import Peppy.Peppy;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.Spectrum;
import Utilities.U;

public class SpectralVisualizer {
	
	public static final int THEME_WHITE = 0;
	public static final int THEME_BLACK = 1;
	private static int theTheme = THEME_WHITE;
	private static boolean cubedScale = true;
	static Color yIonColor = Color.RED;
	static Color bIonColor = Color.BLUE;
	static Color noIonColor = Color.gray;
	static Color bothIonColor = Color.green;
	static Color bkgndColor = Color.white;
	static Color axesColor = Color.black;
	
//	public static void main(String args[]) {
//		if (theTheme == THEME_BLACK) {
//			yIonColor = Color.yellow;
//			bIonColor = Color.cyan;
//			noIonColor = Color.gray;
//			bothIonColor = Color.green;
//			bkgndColor = Color.black;
//			axesColor = Color.white;
//		}
//		if (theTheme == THEME_WHITE) {
//			yIonColor = Color.RED;
//			bIonColor = Color.BLUE;
//			noIonColor = Color.LIGHT_GRAY;
//			bothIonColor = Color.green;
//			bkgndColor = Color.white;
//			axesColor = Color.black;
//		}
//		U.p("drawing images for spectra");
//		Peppy.init();
//		ArrayList<Spectrum> spectra = Spectrum.loadSpectraFromFolder(Properties.spectraDirectoryOrFile);
//		ArrayList<Peptide> peptides = loadPeptides("peptides.txt");
//		generateFullReport(spectra, peptides);
//		U.p("done");
//	}
	
	public static void main(String args[]) {
		System.setProperty("java.awt.headless", "true"); 
		if (theTheme == THEME_BLACK) {
			yIonColor = Color.yellow;
			bIonColor = Color.cyan;
			noIonColor = Color.gray;
			bothIonColor = Color.green;
			bkgndColor = Color.black;
			axesColor = Color.white;
		}
		if (theTheme == THEME_WHITE) {
			yIonColor = Color.RED;
			bIonColor = Color.BLUE;
			noIonColor = Color.LIGHT_GRAY;
			bothIonColor = Color.green;
			bkgndColor = Color.white;
			axesColor = Color.black;
		}
		U.p("drawing images for spectra");
		Peppy.init();
		
		//load spectra
		U.p("loading spectra...");
//		ArrayList<Spectrum> spectra = Spectrum.loadSpectraFromFolder("spectra encode membrane/GO_mem_FASP_dta20100628");
		ArrayList<Spectrum> spectra = Spectrum.loadSpectraFromFolder("spectra encode membrane/SDS");
		
		//load list of spectra and peptides
		File listFolder = new File("differences_spectrum_reports/sds");
		File [] lists = listFolder.listFiles();
		for (int i = 0; i < lists.length; i++) {
			if (!lists[i].getName().endsWith(".txt")) continue;
			U.p("loading list " + lists[i].getName());
			ArrayList<String> spectraNames = new ArrayList<String>();
			ArrayList<String> peptideStrings = new ArrayList<String>();
			ArrayList<String> eValues = new ArrayList<String>();
			try {
				BufferedReader br = new BufferedReader(new FileReader(lists[i]));
				String line = br.readLine();
				int lineNumber = 1;
				while (line != null) {
					String [] chunks = line.split("\t");
					if (chunks == null) {
						U.p("chunks is null");
						U.p("line number: " + lineNumber);
						break;
					}
					if (chunks.length < 3) {
						U.p("chunks length isn't 3");
						U.p("line number: " + lineNumber);
						break;
					}
					spectraNames.add(chunks[0]);
					peptideStrings.add(chunks[1]);
					eValues.add(chunks[2]);
					line = br.readLine();
				}
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
			//find the spectra from the list
			U.p("finding spectra from list...");
			ArrayList<Spectrum> spectraSelect = new ArrayList<Spectrum>();
			ArrayList<Peptide> peptides = new ArrayList<Peptide>();
			for (int j = 0; j < spectraNames.size(); j++) {
				String name = spectraNames.get(j);
				for (Spectrum spectrum: spectra) {
					if (spectrum.getFile().getName().equals(name)) {
						spectraSelect.add(spectrum);
						peptides.add(new Peptide(peptideStrings.get(j)));
						break;
					}
				}
			}
			
			//make the report
			U.p("making drawings...");
			SpectralVisualizer.generateFullReport(spectraSelect, peptides, eValues, lists[i].getName().substring(0, lists[i].getName().length() - 4));
		}

		U.p("done");
	}
	
	public static void generateFullReport(ArrayList<Spectrum> spectra, ArrayList<Peptide> peptides, ArrayList<String> eValues, String folderName) {
		File reportFolder = new File(Properties.reportDirectory, folderName);
		reportFolder.mkdirs();
		File reportFile = new File(reportFolder, "index.html");
		try {
			//create our report directories
			File imageFolder = new File(reportFolder, "images");
			imageFolder.mkdirs();
			
			//set up our main index file
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(reportFile)));
			
			//print headers
			U.appendFile(pw, Properties.reportWebHeaderFile);
			
			for (int i = 0 ; i <spectra.size(); i++) {
				//mark potential peptide peaks
				markMatchingIons(spectra.get(i), peptides.get(i));
				
				//make an image file of the spectrum
				File imageFile = new File(imageFolder, i + ".jpg");
//				drawSpectrum(spectra.get(i), 1000, 300, imageFile);
				drawDeluxSpectrum(spectra.get(i), peptides.get(i), imageFile);
				
				//include in report file
				pw.println("<img src=\"images/" + i + ".jpg\">");
				pw.println("<br>");
				pw.println(peptides.get(i).getAcidSequence());
				pw.println("<br>");
				pw.println(eValues.get(i));
				pw.println("<p>");
			}
			
			//print headers
			U.appendFile(pw, Properties.reportWebHeaderFile);

			pw.flush();
			pw.close();		
		} catch (FileNotFoundException e) {
			U.p("could not find file: " + reportFile.getName());
			e.printStackTrace();
		} catch (IOException e) {
			U.p("could not read file: " + reportFile.getName());
			e.printStackTrace();
		}
	}
	
	public static void generateFullReport(ArrayList<Spectrum> spectra, ArrayList<Peptide> peptides) {
		File reportFolder = new File(Properties.reportDirectory, "spectral visualizations " + System.currentTimeMillis());
		reportFolder.mkdirs();
		File reportFile = new File(reportFolder, "index.html");
		try {
			//create our report directories
			File imageFolder = new File(reportFolder, "images");
			imageFolder.mkdirs();
			
			//set up our main index file
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(reportFile)));
			
			//print headers
			U.appendFile(pw, Properties.reportWebHeaderFile);
			
			for (int i = 0 ; i <spectra.size(); i++) {
				//mark potential peptide peaks
				markMatchingIons(spectra.get(i), peptides.get(i));
				
				//make an image file of the spectrum
				File imageFile = new File(imageFolder, i + ".jpg");
//				drawSpectrum(spectra.get(i), 1000, 300, imageFile);
				drawDeluxSpectrum(spectra.get(i), peptides.get(i), imageFile);
				
				//include in report file
				pw.println("<img src=\"images/" + i + ".jpg\">");
				pw.println("<br>");
				pw.println(peptides.get(i).getAcidSequence());
				pw.println("<p>");
			}
			
			//print headers
			U.appendFile(pw, Properties.reportWebHeaderFile);

			pw.flush();
			pw.close();		
		} catch (FileNotFoundException e) {
			U.p("could not find file: " + reportFile.getName());
			e.printStackTrace();
		} catch (IOException e) {
			U.p("could not read file: " + reportFile.getName());
			e.printStackTrace();
		}
	}
	
	public static void drawDeluxSpectrum(Spectrum spectrum, Peptide peptide, File imageFile)  throws IOException {
		
		//getting sectrum image
		int spectrumWidth = 1000;
		int spectrumHeight = 300;
		BufferedImage spectrumImage = new BufferedImage(spectrumWidth, spectrumWidth, BufferedImage.TYPE_INT_RGB);
		Graphics2D spectrumGraphics = spectrumImage.createGraphics();
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
//		ColorConvertOp cco = new ColorConvertOp(spectrumGraphics.getRenderingHints());
//		chartGraphics.drawImage(chartImage, null, 50, 50);
		
		//fill in a stupid black rectangle
		chartGraphics.fillRect(chartX,chartY + spectrumHeight + 3, spectrumWidth,chartHeigth - chartY - spectrumHeight);
		
		//draw axes
		chartGraphics.setColor(axesColor);
		chartGraphics.drawLine(chartX, chartY + spectrumHeight + 2, chartX + spectrumWidth, chartY + spectrumHeight + 2); //x axis
		chartGraphics.drawLine(chartX + spectrumWidth, chartY + spectrumHeight -10, chartX + spectrumWidth, chartY + spectrumHeight +10); //x axis cap
		chartGraphics.drawLine(chartX, chartY, chartX, chartY + spectrumHeight); //y axis
		chartGraphics.drawLine(chartX - 10, chartY, chartX + 10, chartY); // y cap
		
		//find x tick marks
		double xLabelIncrement = Math.log10(spectrum.getPrecursorMass());
		xLabelIncrement = Math.pow(10, Math.round(xLabelIncrement) - 1);
		int xPixelIncrement = (int) ((xLabelIncrement / spectrum.getPrecursorMass()) * spectrumWidth);
		
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
	
	public static ArrayList<Peptide> loadPeptides(String peptideFileName) {
		ArrayList<Peptide> peptides = new ArrayList<Peptide>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(peptideFileName));
			String line = br.readLine();
			while (line != null) {
				line = line.trim();
				if (!line.equals("")) {
					peptides.add(new Peptide(line));
				}
				line = br.readLine();
			}
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return peptides;
	}
	
	/**
	 * defaults to drawing the ion labels.
	 * @param spectrum
	 * @param width
	 * @param height
	 * @param dest
	 * @throws IOException
	 */
	public static void drawSpectrum(Spectrum spectrum, int width, int height, File dest) throws IOException {
		//setting up Graphics context
		BufferedImage bdest = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		drawSpectrum(spectrum, width, height, dest, true, "", "", bdest);
	}
	
	public static void drawSpectrum(Spectrum spectrum, int width, int height, File dest, boolean drawLabels) throws IOException {
		//setting up Graphics context
		BufferedImage bdest = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		drawSpectrum(spectrum, width, height, dest, drawLabels, "", "", bdest);
	}
	
	public static void drawSpectrum(Spectrum spectrum, int width, int height, File dest, boolean drawLabels, String message1, String message2) throws IOException {
		BufferedImage bdest = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		drawSpectrum(spectrum, width, height, dest, drawLabels, message1, message2, bdest);
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
		double maxValue = spectrum.getPrecursorMass();
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
//					if (peak.getbIonNumber() != -1) {
//						int y1 = height - (height * peak.getbIonNumber() / 20);
//						int labelWidth = fontMetrics.stringWidth(label.toString());
//						int x1 = labelWidth + 4;
//						g.setStroke(new BasicStroke(1.5f));
//						g.drawLine(x1, y1, xLoc, yLoc);
//						g.drawString("" + label, 2, y1);
//					}
//					if (peak.getyIonNumber() != -1) {
//						int y1 = height - (height * peak.getyIonNumber() / 20);
//						int labelWidth = fontMetrics.stringWidth(label.toString());
//						int x1 = width - (labelWidth + 4);
//						g.setStroke(new BasicStroke(1.5f));
//						g.drawLine(x1, y1, xLoc, yLoc);
//						g.drawString("" + label, width - labelWidth, y1);
//					}
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
	

	
	public static void markMatchingIons(Spectrum spectrum, Peptide peptide) {
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
		theoreticalPeakMass = peptide.getMass() + Properties.rightIonDifference;
		for (i = 0; i < peptideLengthMinusOne; i++) {
			theoreticalPeakMass -= AminoAcids.getWeightMono(acidSequence[i]);
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
					spectrum.getPeak(peakIndex).setColor(yIonColor);
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
					if (spectrum.getPeak(peakIndex).getColor().equals(yIonColor)) {
						spectrum.getPeak(peakIndex).setColor(bothIonColor);
					} else {
						spectrum.getPeak(peakIndex).setColor(bIonColor);
					}
					spectrum.getPeak(peakIndex).setbIonNumber(seqIndex + 1);
				}
			}
			
			peakIndex++;
		}
		
	}
	

}
