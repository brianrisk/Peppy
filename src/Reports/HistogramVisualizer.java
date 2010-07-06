package Reports;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

public class HistogramVisualizer {
	
	/**
	 * makes a histogram which has as many pixels for width the histogram array length
	 * 
	 * @param histogram
	 * @param height
	 * @param dest
	 * @throws IOException
	 */
	public static void drawHistogram(int [] histogram, int height, File dest) throws IOException {
		
		int width = histogram.length;
		
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
		int top;
		double scaleFactor = (double) height / max;
		for (int i = 0; i < histogram.length; i++) {
			top = height - (int) (histogram[i] * scaleFactor);
			g.drawLine(i, top, i, height);
			
		}
		ImageIO.write(bdest,"JPG",dest);
	}

}
