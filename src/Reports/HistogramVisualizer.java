package Reports;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

import Utilities.U;

public class HistogramVisualizer {
	
	/**
	 * makes a histogram which has as many pixels for width the histogram array length
	 * 
	 * @param histogram
	 * @param height
	 * @param dest
	 * @throws IOException
	 */
	public static void drawShadedHistogram(int [] histogram, int height, int maxBar, File dest) throws IOException {
		
		int width = histogram.length;
		
		//setting up Graphics context
		BufferedImage bdest = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		Graphics2D g = bdest.createGraphics();
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_SPEED));
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF));
		g.setColor(Color.white);
		g.fillRect(0,0,width,height);
		
		//find the maximum value
//		int max = 0;
//		for (int i = 0; i < histogram.length; i++) {
//			if (histogram[i] > max) max = histogram[i];
//		}
		
		//draw the bars
		int top;
		double scaleFactor = (double) 255 / 100;
		for (int i = 0; i < histogram.length; i++) {
			
			top = (int) (histogram[i] * scaleFactor);
			if (top > 0){
				if (top > 255) top = 255;
				g.setColor(new Color(0, 0, 0, top));
				g.drawLine(i, 0, i, height);
			}
			
		}
		ImageIO.write(bdest,"JPG",dest);
	}
	
	public static void drawHistogram(int [] histogram, int height, int width, File dest) throws IOException {
		//setting up Graphics context
		BufferedImage bdest = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		Graphics2D g = bdest.createGraphics();
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_SPEED));
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF));
		g.setColor(Color.white);
		g.fillRect(0,0,width,height);
		
		//find the maximum value
		int max = 0;
		for (int i = 0; i < histogram.length; i++) {
			if (histogram[i] > max) max = histogram[i];
		}
		
		//draw the bars
		g.setColor(Color.darkGray);
		g.drawRect(0,0,width - 1, height - 1);
		double scaleFactor = (double) height / max;
		for (int i = 0; i < histogram.length; i++) {
			int x1 = (int) ((double) i * width / histogram.length);
			int barWidth = (int) ((double) (i + 1.0) * width / histogram.length) - x1;
			if (barWidth < 1) barWidth = 1;
			int barHeight = (int) (histogram[i] * scaleFactor);
			//account for tiny bars that might not be displayed
			if (barHeight <= 1 && histogram[i] > 0) {
				barHeight = 5;
				g.setColor(Color.lightGray);
			} else {
				g.setColor(Color.black);
			}
			g.fillRect(x1, height - barHeight, barWidth, barHeight);
		}
		
		//outline
		g.setColor(Color.black);
		g.drawRect(0,0,width - 1, height - 1);
		
		/*
		 * calculate and draw normal curve
		 */
		//find  mean
		double variance = 0;
		double mean = 0;
		int total = 0;
		for (int i = 0; i < histogram.length; i++) {
			mean += i * histogram[i];
			total += histogram[i];
		}
		mean /= total;
		//draw mean line
		g.setColor(Color.green);
		int meanX = (int) ((double) mean * width / histogram.length);
		g.drawLine(meanX, 0, meanX, height);
		//find variance
		double meanDifference;
		for (int i = 0; i < histogram.length; i++) {
			meanDifference = i - mean;
			variance += histogram[i] * (meanDifference * meanDifference);
		}
		variance /= total;
		
		//draw normal curve
		double normalConstant = 1.0 / Math.sqrt(2 * Math.PI * variance);
		double numerator;
		double denominator = -2 * variance;
		g.setColor(Color.red);
		g.setStroke(new BasicStroke(2.0f));
		int x1 = 0, y1 = height, x2, y2;
		
		
		for (int i = 0; i < histogram.length; i++) {
			numerator = (i - mean);
			numerator *= numerator;
			x2 = (int) ((double) i * width / histogram.length);
			y2 = height - (int) (scaleFactor * total * normalConstant * Math.exp(numerator / denominator));
			g.drawLine(x1, y1, x2, y2);
			x1 = x2;
			y1 = y2;
		}
		
		//draw log normal curve
		g.setColor(Color.green);
		g.setStroke(new BasicStroke(2.0f));
		x1 = 0;
		y1 = height;
		double logNormalFirstTerm;
//		if (
//				dest.getName().startsWith("107") ||
//				dest.getName().startsWith("107") ||
//				dest.getName().startsWith("107") ||
//				dest.getName().startsWith("107") 
//		) {
//			U.p("variance: " + variance);
//		}
		for (int i = 0; i < histogram.length; i++) {
			logNormalFirstTerm = 1.0 / (i * Math.sqrt(2 * Math.PI * variance));
			numerator = Math.log(i) - mean;
			numerator *= numerator;
			x2 = (int) ((double) i * width / histogram.length);
			y2 = height - (int) (scaleFactor * total * logNormalFirstTerm * Math.exp(numerator / denominator));
			g.drawLine(x1, y1, x2, y2);
			x1 = x2;
			y1 = y2;
		}
		
		
		
		
		
		//write
		ImageIO.write(bdest,"JPG",dest);
	}
	
	public static void drawHistogram(double [] histogram, int height, int width, File dest) throws IOException {
		//setting up Graphics context
		BufferedImage bdest = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		Graphics2D g = bdest.createGraphics();
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_SPEED));
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF));
		g.setColor(Color.white);
		g.fillRect(0,0,width,height);
		
		//find the maximum value
		double max = 0;
		for (int i = 0; i < histogram.length; i++) {
			if (histogram[i] > max) max = histogram[i];
		}
		
		//draw the bars
		g.setColor(Color.black);
		g.drawRect(0,0,width - 1, height - 1);
		double scaleFactor = (double) height / max;
		boolean redLine = false;
		for (int i = 0; i < histogram.length; i++) {
			int x1 = (int) ((double) i * width / histogram.length);
			int barWidth = (int) ((double) (i + 1.0) * width / histogram.length) - x1;
			if (barWidth < 1) barWidth = 1;
			int barHeight = (int) (histogram[i] * scaleFactor);
			
//			U.p(x1 + ", " + (height - barHeight)+ ", " +  barWidth+ ", " +  barHeight);
			g.fillRect(x1, height - barHeight, barWidth, barHeight);
			if (histogram[i] <= 0.1 && redLine == false) {
				redLine = true;
				g.setColor(Color.red);
				g.drawLine(x1, 0, x1, height);
			}
		}
		ImageIO.write(bdest,"JPG",dest);
	}

}
