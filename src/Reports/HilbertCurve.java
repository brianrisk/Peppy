package Reports;

import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;

import Peppy.Match;
import Utilities.U;

public class HilbertCurve {
	
	public static final int UP = 0;
	public static final int DOWN = 1;
	public static final int LEFT = 2;
	public static final int RIGHT = 3;

	private int curveLength;
	private int [] powers;
	private int sideLength;
	private int exponent;
	private BufferedImage curveImage;
	
	private final int [][][] vectors = {
			{{0,0,1,1},
			 {0,1,1,0}},
			{{1,1,0,0},
			 {1,0,0,1}},
			{{1,0,0,1},
			 {1,1,0,0}},
			{{0,1,1,0},
			 {0,0,1,1}}
	};
	
	private static final int [][] directions = {
			{RIGHT, UP, UP, LEFT},
			{LEFT, DOWN, DOWN, RIGHT},
			{UP, LEFT, LEFT, DOWN},
			{DOWN, RIGHT, RIGHT, UP}	
	};
	
	public static void main(String args[]) {
		new HilbertCurve(0, 1000);
	}
	
	/**
	 * 
	 * @param match
	 * @param sequenceIndexStart inclusive
	 * @param sequenceIndexStop exclusive. that is, if value is 1000 then the max index value is 999.
	 */
	public HilbertCurve(int sequenceIndexStart, int sequenceIndexStop) {
		//find the length of the side of the Hilbert Curve
		double workLength = Math.sqrt(sequenceIndexStop - sequenceIndexStart);
		workLength = U.log(2, workLength);
		workLength = Math.ceil(workLength);
		exponent = (int) workLength;
		sideLength = (int) Math.pow(2, exponent);
		curveLength = sideLength * sideLength;
		
		//store so we don't have to repeatedly do these divisions by two
		powers = new int[exponent];
		int power = sideLength / 2;
		for (int i = 0; i < exponent; i++) {
			powers[i] = power;
			U.p(power);
			power /= 2;
		}
		Point point = getPoint(1);
		U.p(point);
		
		//initialize our image
//		curveImage = new BufferedImage(sideLength, sideLength, BufferedImage.TYPE_INT_RGB);
	}
	
	public Point getPoint(int index) {
		int x = 0;
		int y = 0;
		int lowerBound = 0;
		int upperBound = curveLength;
		int quadrant;
		int myDirection = UP;
		
		int quadrantSize = curveLength / 4;
		for (int i = 0; i < exponent; i++) {
			U.p("ub:" + upperBound + ", lb:" + lowerBound);
			if (lowerBound == index) break;
			quadrant = getQuadrant(index, lowerBound, upperBound);
			U.p("quadrant: " + quadrant);
			x += vectors[myDirection][0][quadrant] * powers[i];
			y += vectors[myDirection][1][quadrant] * powers[i];
			lowerBound += quadrant * quadrantSize;
			upperBound -= (3 - quadrant) * quadrantSize;
			quadrantSize /= 4;
			myDirection = directions[myDirection][quadrant];
		}
		return new Point(x, y);
	}
	
	private int getQuadrant(int index, int lowerBound, int upperBound) {
		int length = (upperBound - lowerBound) / 4;
		int boundary = lowerBound + length;
		for (int i = 0; i < 4; i++) {
			if (index < boundary) return i;
			boundary += length;
		}
		return -1;
	}
	


	



}
