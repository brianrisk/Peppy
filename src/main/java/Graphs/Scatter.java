package Graphs;

import Peppy.U;

import java.util.ArrayList;

/**
 * Contains static methods for calculating least squares
 * @author Brian Risk
 *
 */
public class Scatter {
	ArrayList<Point> points;
	
	private double slope;
	private double intercept;
	private double correlation;
	private boolean leastSquareHasBeenCalculated = false;
	
	public Scatter() {
		points = new ArrayList<Point>();
	}
	
	public Scatter(ArrayList<Point> points) {
		this.points = points;
	}
	
	public void addPoint(Point point) {
		points.add(point);
	}
	
	/* constructor used for e-value, p-value */
	public Scatter(double [] xValues, double [] yValues, int [] histogram, int start, int stop) {
		points = new ArrayList<Point>();
		for (int i = start; i < stop; i++) {
			if (histogram[i] > 0) {
				points.add(new Point(xValues[i], yValues[i]));
			}
		}
	}
	

	public void calculateLeastSquare() {
		double meanX = 0;
		double meanY = 0;
		double meanXY = 0;
		double meanXX = 0;
		double meanYY = 0;
		for (Point point: points) {
			meanX += point.getX();
			meanY += point.getY();
			meanXY += point.getX() * point.getY();
			meanXX += point.getX() * point.getX();
			meanYY += point.getY() * point.getY();
		}
		meanX /= points.size();
		meanY /= points.size();
		meanXY /= points.size();
		meanXX /= points.size();
		meanYY /= points.size();
		
		slope = (meanXY - meanX * meanY) / (meanXX - meanX * meanX);
		intercept = meanY - slope * meanX;
		correlation = slope * (meanXY - meanX * meanY) / (meanYY - meanY * meanY);	
		
		/* mark that this has been calculated */
		leastSquareHasBeenCalculated = true;
	}
	
	public double calculatePerfectP() {
		/* find boundaries */
		double maxX = Double.NEGATIVE_INFINITY;
		double minX = Double.MAX_VALUE;
		double maxY = Double.NEGATIVE_INFINITY;
		double minY = Double.MAX_VALUE;
		for (Point point: points){
			if (point.getX() > maxX) maxX = point.getX();
			if (point.getX() < minX) minX = point.getX();
			
			if (point.getY() > maxY) maxY = point.getY();
			if (point.getY() < minY) minY = point.getY();
		}
		
		/* ranges */
		double xRange = maxX - minX;
		double yRange = maxY - minY;
		
		/* get the combined multiples of the probabilities */
		double p = 1;
		double area;
		for (Point point: points) {
			area = 2 * Math.abs(point.x - point.y);
			U.p(area);
			p *= (area / yRange);
		}
		
		return p;
	}
	

	public double getSlope() {
		if (!leastSquareHasBeenCalculated) calculateLeastSquare();
		return slope;
	}

	public double getIntercept() {
		if (!leastSquareHasBeenCalculated) calculateLeastSquare();
		return intercept;
	}

	public double getCorrelation() {
		if (!leastSquareHasBeenCalculated) calculateLeastSquare();
		return correlation;
	}
	
	public ArrayList<Point> getPoints() {
		return points;
	}
	

}
