package Math;

import java.awt.geom.Point2D;
import java.util.ArrayList;

/**
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class LeastSquares {

	public static double calculateM( double [] yValues, int [] histogram, int start, int stop) {
		double numerator1, numerator2, denomenator1, denomenator2;
		double numerator, denomenator;
		double temp1 = 0.0, temp2 = 0.0, temp = 0.0;
		double parameterM;
		int i;
		for (i = start; i < stop; i++) {
			if (histogram[i] > 0)
				temp += (i * yValues[i]);
		}
		numerator1 = (stop - start) * (temp);
		for (i = start; i < stop; i++) {
			if (histogram[i] > 0)
				temp1 += i;
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
				temp1 = i;
				temp2 += (temp1 * temp1);
			}
		}
		denomenator1 = (stop - start) * temp2;

		temp1 = 0.0; 
		for (i = start; i < stop; i++) {
			if (histogram[i] > 0)
				temp1 += i;
		}
		denomenator2 = (temp1 * temp1);
		denomenator = denomenator1 - denomenator2;
		parameterM = numerator / denomenator;
		return parameterM;
	}
	
	public static double calculateB( double [] yValues, int [] histogram, int start, int stop, double m) {
		double parameterB;
		double temp1, temp2;
		int i;

		temp1 = temp2 = 0.0;
		for (i = start; i < stop; i++) {
			if (histogram[i] > 0)
				temp1 += i;
		}

		for (i = start; i < stop; i++) {
			if (histogram[i] > 0)
				temp2 += yValues[i];
		}
		parameterB = (1.0 / (stop - start)) * (temp2 - m * temp1);
		return parameterB;
	}


}
