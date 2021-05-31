package Experimental;

import Peppy.AminoAcids;
import Peppy.Peppy;
import Peppy.Properties;
import Peppy.U;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.*;

public class DrawProbabilityGrid {
	
	public static void main(String args[]) {
		/* set up initial state */
		Peppy.init(args);
		
		int rows = (int) (Properties.peptideMassMaximum - Properties.peptideMassMinimum + 1);
		int cols = (int) ((Properties.peptideMassMaximum - Properties.peptideMassMinimum) / AminoAcids.getWeightMono(AminoAcids.G) );
		double [][] results = new double[rows][cols];
		
		/* load in the results matrix */
		try {
			BufferedReader br = new BufferedReader(new FileReader("CommonIonCount Human.txt"));
			String line = br.readLine();
			int rowIndex = 0;
			while (line != null) {
				String [] chunks = line.split("\t");
				for (int colIndex = 0; colIndex < chunks.length; colIndex++) {
					results[rowIndex][colIndex] = Integer.parseInt(chunks[colIndex]);
				}
				
				line = br.readLine();
				rowIndex++;
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		/* take the log! */
		double maxLog = 0;
		for (int rowIndex = 0; rowIndex < rows; rowIndex++) {
			for (int colIndex = 0; colIndex < cols; colIndex++) {
				if (results[rowIndex][colIndex] != 0) {
					results[rowIndex][colIndex] = Math.log(results[rowIndex][colIndex]);
					results[rowIndex][colIndex] += 1;
					if (results[rowIndex][colIndex] > maxLog) maxLog = results[rowIndex][colIndex];
				}
			}
		}
		
		/* create percent matrix */
		double [][] percentages = new double[rows][cols];
		
		
		for (int rowIndex = 0; rowIndex < rows; rowIndex++) {
			double rowSum = 0;
			double rowMax = 0;
			for (int colIndex = 0; colIndex < cols; colIndex++) {
				rowSum += results[rowIndex][colIndex];
				if (results[rowIndex][colIndex] > rowMax) rowMax = results[rowIndex][colIndex];
			}
			double normalizingFactor = (double) rowSum / rowMax;
			for (int colIndex = 0; colIndex < cols; colIndex++) {
				double percentage = 0;
				if (rowSum != 0) {
					percentage = (double) results[rowIndex][colIndex] / maxLog;
//					percentage *= normalizingFactor;
					if (percentage > 1) percentage = 1;
				}
				percentages[rowIndex][colIndex] = percentage;
			}
		}
		
		
		/* draw the grid */
		int cellWidth = 10;
		int cellHeight = 1;
		int verticalScale = 4;
		cols = 90;
//		rows = 000;
		int imageWidth = cellWidth * cols;
		int imageHeight = cellHeight * rows / verticalScale;
		BufferedImage bdest = new BufferedImage(imageWidth, imageHeight, BufferedImage.TYPE_INT_RGB);
		Graphics2D g = bdest.createGraphics();
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_SPEED));
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF));
		for (int rowIndex = 0; rowIndex < rows; rowIndex += verticalScale) {
			for (int colIndex = 0; colIndex < cols; colIndex++) {
				Color color;
				if (percentages[rowIndex][colIndex] == 0) {
					color = Color.black;
				} else {
					color = Color.getHSBColor((float) percentages[rowIndex][colIndex], 1, 1);
				}
				int x = colIndex * cellWidth;
				int y = rowIndex * cellHeight / verticalScale;
				g.setColor(color);
				g.fillRect(x, y, cellWidth, cellHeight);
			}
		}
		try {
			ImageIO.write(bdest,"PNG",new File("probailityGrid.png"));
		} catch (IOException e) {
			e.printStackTrace();
		}
		U.p("done");
	}
	

}
