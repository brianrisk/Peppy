package Graphs;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

import javax.imageio.ImageIO;

import Peppy.U;

public class ScatterVisualizer {
	
	Scatter scatter;
	File file;
	private static int width = 1000;
	

	double max = Double.MIN_VALUE;
	double min;
	
	/* test it out! */
	public static void main(String[] args) {
		File file = new File("scatter test.png");
		Random random = new Random();
		ArrayList<Point> points = new ArrayList<Point>();
		for (int i = 0; i < 100000; i++) {
			points.add(new Point(random.nextGaussian(), random.nextGaussian()));
		}
		Scatter scatter = new Scatter(points);
		ScatterVisualizer me = new ScatterVisualizer(scatter, file);
		me.draw();
		
		U.p("slope: " + scatter.getSlope());
		U.p("intercept: " + scatter.getIntercept());
		U.p("correlation: " + scatter.getCorrelation());
//		U.p("perfect correlation P value: " + scatter.calculatePerfectP());
		
		U.p("done");
	}
	
	
	public ScatterVisualizer(Scatter scatter, File file) {
		super();
		this.scatter = scatter;
		this.file = file;
		
		/* find boundaries */
		double value;
		for (Point point: scatter.getPoints()){
			value = Math.abs(point.x);
			if (value > max) max = value;
			value = Math.abs(point.y);
			if (value > max) max = value;
		}
		min = -max;
	}
	
	public void draw() {
		/* setting up Graphics context */
		BufferedImage bufferedImage = new BufferedImage(width, width, BufferedImage.TYPE_INT_RGB);
		Graphics2D g = bufferedImage.createGraphics();
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY));
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
		g.setColor(Color.white);
		g.fillRect(0,0,width,width);
		
		/* draw the axes */
		g.setColor(Color.black);
		int yZero = getScaledValue(0,true);
		g.drawLine(0, yZero, width, yZero);
		int xZero = getScaledValue(0,false);
		g.drawLine(xZero, 0, xZero, width);
		
		/* draw the points */
		g.setColor(new Color(0,0,0, 255));
		g.setColor(new Color(0,0,128, 255));
		for (Point point: scatter.getPoints()) {
//			g.drawLine(
//					getScaledValue(point.getX(), false),
//					getScaledValue(point.getY(), true),
//					getScaledValue(point.getX(), false),
//					getScaledValue(point.getY(), true)
//					);
			g.drawOval(
					getScaledValue(point.getX(), false) - 2,
					getScaledValue(point.getY(), true) - 2,
					5,
					5
					);
		}
		
		/* draw the best fit slope */
		g.setColor(new Color(200,0,0, 128));
		g.setStroke(new BasicStroke(2.0f));
		double y1 = scatter.getSlope() * min + scatter.getIntercept();
		double y2 = scatter.getSlope() * max + scatter.getIntercept();
		g.drawLine(0, getScaledValue(y1, true), width, getScaledValue(y2, true));
		
		
		/* draw perfect correlation slope */
		g.setColor(new Color(0,0,200, 128));
		g.setStroke(new BasicStroke(2.0f));
		g.drawLine(0, getScaledValue(min, true), width, getScaledValue(max, true));
		
		
		/* write */
		try {
			ImageIO.write(bufferedImage,"PNG",file);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private int getScaledValue(double value, boolean isY) {
		double slope = (width) / (max - min);
		double intercept = min * slope * -1;
		double location = value * slope + intercept;
		if (isY) location = width - location;
		return (int) location;
	}
		

	


}
