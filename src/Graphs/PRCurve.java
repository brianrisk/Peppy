package Graphs;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.Polygon;
import java.awt.RenderingHints;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;

import javax.imageio.ImageIO;

public class PRCurve {
	
	private double areaUnderCurve = 0;
	ArrayList<Point2D.Double> points;
	int width = 500;
	Color color = Color.BLACK;
	
	
	public PRCurve(ArrayList<Point2D.Double> points) {
		this.points = points;
		
		/* adjust for non-monotonic decreasing PR curve */
		double previousPrecision = points.get(points.size() - 1).x;
		for (int i = points.size() - 2; i >= 0; i--) {
			Point2D.Double point = points.get(i);
			if (point.x < previousPrecision) {
				point.x = previousPrecision;
			} else {
				previousPrecision = point.x;
			}
		}
		calculateAreaUnderCurve();
	}
	
	public double calculateAreaUnderCurve() {
		areaUnderCurve = 0;
		double recall;
		double recallPrevious = points.get(0).y;
		for (int i = 1; i < points.size(); i++) {
			recall =  points.get(i).y;
			areaUnderCurve += (recall - recallPrevious) *  points.get(i).x;
			recallPrevious = recall;
		}
		return areaUnderCurve;
	}
	
	public BufferedImage generatePrecisionRecallCurve() {

		/* setting up Graphics context */
		BufferedImage bufferedImage = new BufferedImage(width, width, BufferedImage.TYPE_INT_RGB);
		Graphics2D g = bufferedImage.createGraphics();
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY));
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
		g.setColor(Color.white);
		g.fillRect(0,0,width,width);
		
		/* construct the polygon */
		Polygon polygon = new Polygon();
		for (Point2D.Double point: points) {
			int x = (int) (point.y * width);
			int y = (int) (width - (point.x * width));
//			U.p(x + ", " + y);
			polygon.addPoint(x, y);
		}
		/* add the corner points */
		polygon.addPoint((int) (points.get(points.size() - 1).y * width), width);
		polygon.addPoint(0, width);
		
		
		/* draw the curve */
		g.setColor(color);
		g.fillPolygon(polygon);
		
		
		//create axis
		int indent = 40;
		int axisSpace = 5;
		BufferedImage axisImage = new BufferedImage(width + indent, width + indent, BufferedImage.TYPE_INT_RGB);
		Graphics2D axisGraphics = axisImage.createGraphics();
		axisGraphics.setFont(new Font("Monospaced", Font.PLAIN, 36));
		FontMetrics fm = axisGraphics.getFontMetrics();
		axisGraphics.setColor(Color.white);
		axisGraphics.fillRect(0, 0, width + indent, width + indent);
		axisGraphics.drawImage(bufferedImage, indent, 0, width, width, null);
		
		//draw axis
		axisGraphics.setColor(Color.black);
		axisGraphics.setStroke(new BasicStroke(3.0f));
		axisGraphics.drawLine(indent, 0, indent, width);
		axisGraphics.drawLine(indent, width, width + indent, width);
		
		//draw axis caps
		int capRadius = 15;
		axisGraphics.drawLine(indent - capRadius, 0, indent + capRadius, 0);
		axisGraphics.drawLine(indent + width, width - capRadius, indent + width, width + capRadius);
		
		//draw labels
		axisGraphics.drawString("0", indent - fm.stringWidth("0") - axisSpace, width + fm.getHeight());
		axisGraphics.drawString("1", indent - fm.stringWidth("1") - axisSpace, fm.getHeight() + axisSpace);
		axisGraphics.drawString("1", width + indent - fm.stringWidth("1") - axisSpace, width + fm.getHeight());
		
		/* label the curve with the percentage of area covered */
		NumberFormat nfPercent = NumberFormat.getPercentInstance();
		nfPercent.setMaximumFractionDigits(2);
		String areaLabel = nfPercent.format(areaUnderCurve);
		axisGraphics.setColor(Color.white);
		axisGraphics.drawString(areaLabel,  indent +  (width - fm.stringWidth(areaLabel)) / 2 , indent +  (width - fm.getHeight()) / 2);
		
		return axisImage;
		
	}

	public void writeFile(File file)  {
		BufferedImage bufferedImage = generatePrecisionRecallCurve();
		try {
			ImageIO.write(bufferedImage,"JPG",file);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
