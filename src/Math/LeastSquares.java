package Math;

import java.awt.geom.Point2D;
import java.util.ArrayList;

/**
 * Contains static methods for calculating least squares
 * @author Brian Risk
 *
 */
public class LeastSquares {
	
	public static double calculateM(double [] xValues, double [] yValues, int start, int stop) {
		return calculateM(xValues, yValues, null, start, stop);
	}
	
	public static double calculateB(double [] xValues, double [] yValues, int start, int stop, double m) {
		return calculateB(xValues, yValues, null, start, stop, m);
	}
	
	public static double calculateM(double [] xValues, double [] yValues, int [] histogram, int start, int stop) {
		double numerator1, numerator2, denomenator1, denomenator2;
		double numerator = 0.0, denomenator = 0.0;
		double temp1 = 0.0, temp2 = 0.0, temp = 0.0;
		double parameterM;
		int i;
		for (i = start; i < stop; i++) {
			if (histogram[i] > 0)
				temp += (xValues[i] * yValues[i]);
		}
		numerator1 = (stop - start) * (temp);
		for (i = start; i < stop; i++) {
			if (histogram[i] > 0)
				temp1 += xValues[i];
		}
		for (i = start; i < stop; i++) {
			if (histogram[i] > 0)
				temp2 += yValues[i];
		}
		numerator2 = temp1 * temp2;
		numerator = numerator1 - numerator2;
		temp1 = temp2 = 0.0;
		for (i = start; i < stop; i++) {
			if (histogram[i] > 0) {
				temp1 = xValues[i];
				temp2 += (temp1 * temp1);
			}
		}
		denomenator1 = (stop - start) * temp2;

		temp1 = 0.0; 
		for (i = start; i < stop; i++) {
			if (histogram[i] > 0)
				temp1 += xValues[i];
		}
		denomenator2 = (temp1 * temp1);
		denomenator = denomenator1 - denomenator2;
		parameterM = numerator / denomenator;
		return parameterM;
	}
	
	public static double calculateB(double [] xValues, double [] yValues, int [] histogram, int start, int stop, double m) {
		double parameterB;
		double temp1, temp2;
		int i;

		temp1 = temp2 = 0.0;
		for (i = start; i < stop; i++) {
			if (histogram[i] > 0)
				temp1 += xValues[i];
		}

		for (i = start; i < stop; i++) {
			if (histogram[i] > 0)
				temp2 += yValues[i];
		}
		parameterB = (1.0 / (stop - start)) * (temp2 - m * temp1);
		return parameterB;
	}
	
	public static double [] calculateLeastSquare(ArrayList<Point2D.Double> points) {
		double meanX = 0;
		double meanY = 0;
		double meanXY = 0;
		double meanXX = 0;
		for (Point2D.Double point: points) {
			meanX += point.getX();
			meanY += point.getY();
			meanXY += point.getX() * point.getY();
			meanXX += point.getX() * point.getX();
		}
		meanX /= points.size();
		meanY /= points.size();
		meanXY /= points.size();
		meanXX /= points.size();
		
		double beta = (meanXY - meanX * meanY) / (meanXX - meanX * meanX);
		
		double alpha = meanY - beta * meanX;
		
		double [] out = {alpha, beta};
		return out;
	}

}
