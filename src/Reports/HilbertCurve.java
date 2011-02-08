package Reports;

import java.awt.Color;
import java.awt.Point;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import javax.imageio.ImageIO;

import Peppy.DNA_DigestionThread;
import Peppy.DNA_Sequence;
import Peppy.Definitions;
import Peppy.Properties;
import Peppy.Sequence;
import Utilities.U;

public class HilbertCurve {
	
	public static final int UP = 0;
	public static final int DOWN = 1;
	public static final int LEFT = 2;
	public static final int RIGHT = 3;

	private int curveLength;
	private int [] quadrantLengths;
	private int [] segmentLengths;
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
			{DOWN, LEFT, LEFT, UP},
			{UP, RIGHT, RIGHT, DOWN}	
	};
	
	/**
	 * draw all of the STOP codons
	 * @param args
	 */
	public static void main(String args[]) {
		U.p("whaa");
		new HilbertCurve();
	}
	
	/**
	 * 
	 * @param match
	 * @param sequenceIndexStart inclusive
	 * @param sequenceIndexStop exclusive. that is, if value is 1000 then the max index value is 999.
	 */
	public HilbertCurve() {
		//get our sequence
		ArrayList<Sequence> sequences = Sequence.loadSequences(Properties.sequenceDirectoryOrFile);
		Sequence sequence = sequences.get(0);
		DNA_Sequence dnaSequence = sequence.getNucleotideSequences().get(0);
		String dna = dnaSequence.getSequence();		
		int sequenceIndexStart = 0;
		int sequenceIndexStop = dna.length();
		
		//find the length of the side of the Hilbert Curve square
		double workLength = Math.sqrt(sequenceIndexStop - sequenceIndexStart);
		workLength = U.log(2, workLength);
		workLength = Math.ceil(workLength);
		exponent = (int) workLength;
		sideLength = (int) Math.pow(2, exponent);
		
		//the length of the entire, unfolded line
		curveLength = sideLength * sideLength;
		
		//store so we don't have to repeatedly do these divisions by two
		quadrantLengths = new int[exponent];
		segmentLengths = new int[exponent];
		int quadrantLength = sideLength / 2;
		int segmentLength = curveLength / 4;
		for (int i = 0; i < exponent; i++) {
			quadrantLengths[i] = quadrantLength;
			segmentLengths[i] = segmentLength;
			quadrantLength /= 2;
			segmentLength /= 4;
		}

		//initialize our image
		int whiteRGB = Color.white.getRGB();
		curveImage = new BufferedImage(sideLength, sideLength, BufferedImage.TYPE_INT_RGB);
		for (int x = 0; x < sideLength; x++) {
			for (int y = 0; y < sideLength; y++) {
				curveImage.setRGB(x, y, whiteRGB);
			}
		}
		
		//draw all STOP points
		char [] codon = new char[3];
		char aminoAcid;
		int mod;
		boolean forwards = true;
		Point point;
		boolean inORF = false;
		Color color;
		for (int start = 0; start < 3; start++) {
			for (int index = start; index < dna.length(); index++) {
				mod = (index - start) % 3;
				codon[mod] = dna.charAt(index);
				if (mod == 2) {
					aminoAcid = Definitions.aminoAcidList[DNA_DigestionThread.indexForCodonArray(codon, forwards)];
					if (aminoAcid == 'M') {
						inORF = true;
					}
					if (aminoAcid == '.') {
						inORF = false;
					}
					
				}
				if (inORF) {
					point = getPoint(index);
					color = new Color(curveImage.getRGB(point.x, point.y));
					color = color.darker();
					curveImage.setRGB(point.x, point.y, color.getRed());
				}
			}
		}
		
		//write the file
		File outFile = new File("ecoli-ORF0rgb.png");
		try {
			ImageIO.write(curveImage,"PNG",outFile);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		U.p("done");
	}
	
//	public static void main(String args[]) {
//		new HilbertCurve(0, 1024);
//	}
//	
//	/**
//	 * 
//	 * @param match
//	 * @param sequenceIndexStart inclusive
//	 * @param sequenceIndexStop exclusive. that is, if value is 1000 then the max index value is 999.
//	 */
//	public HilbertCurve(int sequenceIndexStart, int sequenceIndexStop) {
//		//find the length of the side of the Hilbert Curve square
//		double workLength = Math.sqrt(sequenceIndexStop - sequenceIndexStart);
//		workLength = U.log(2, workLength);
//		workLength = Math.ceil(workLength);
//		exponent = (int) workLength;
//		sideLength = (int) Math.pow(2, exponent);
//		
//		//the length of the entire, unfolded line
//		curveLength = sideLength * sideLength;
//		
//		//store so we don't have to repeatedly do these divisions by two
//		quadrantLengths = new int[exponent];
//		segmentLengths = new int[exponent];
//		int quadrantLength = sideLength / 2;
//		int segmentLength = curveLength / 4;
//		for (int i = 0; i < exponent; i++) {
//			quadrantLengths[i] = quadrantLength;
//			segmentLengths[i] = segmentLength;
//			quadrantLength /= 2;
//			segmentLength /= 4;
//		}
//		
//		//initialize our image
//		curveImage = new BufferedImage(sideLength, sideLength, BufferedImage.TYPE_INT_RGB);
//		
//		for (int i = 0; i < 1024; i++) {
//			Point point = getPoint(i);
//			curveImage.setRGB(point.x, sideLength - point.y - 1, Color.WHITE.getRGB());
//			File imageFile = new File("hilbert/" + i + ".png");
//			try {
//				ImageIO.write(curveImage,"PNG",imageFile);
//			} catch (IOException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//		}
//		U.p("done");
//	}
	
	public Point getPoint(int index) {
		int x = 0;
		int y = 0;
		int lowerBound = 0;
		int upperBound = curveLength;
		int quadrant = 0;
		int myDirection = UP;
		for (int i = 0; i < exponent; i++) {
			quadrant = getQuadrant(index, lowerBound, upperBound);
			x += vectors[myDirection][0][quadrant] * quadrantLengths[i];
			y += vectors[myDirection][1][quadrant] * quadrantLengths[i];
			lowerBound += quadrant * segmentLengths[i];
			upperBound -= (3 - quadrant) * segmentLengths[i];
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
