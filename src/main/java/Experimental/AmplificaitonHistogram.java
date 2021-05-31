package Experimental;

import Navigator.Match;
import Navigator.Report;
import Peppy.U;

import java.io.*;
import java.util.ArrayList;

/**
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class AmplificaitonHistogram extends Report {
	
	


	private static int bucketSize = 2000000;
	private static double maxBarValue = 5;
	
	public AmplificaitonHistogram(ArrayList<Match> matches) {
		super(matches);
	}
	
	
	
	private static void addMatchesToHistogram(ArrayList<Match> matches, double [] histogram, String sequenceName) {
		int bucketIndex;
		double increment = 100000.0 / matches.size();
		for (Match match: matches) {
			if (match.getString("sequenceName").equals(sequenceName)) {
				int matchCount =  match.getInt("RankCount");
				int start = match.getInt("start");
				bucketIndex = start / bucketSize;
				histogram[bucketIndex] += increment / matchCount;
			}
		}
	}
	
	
	public File createReport(File reportDirectory) {
		
		File outputFile = new File(reportDirectory, "expression histogram.txt");
		
		/* set up histogram */
		String sequenceName = "chr8";
		int chrLength = 146364022;
		int numberOfBuckets = chrLength / bucketSize + 1;
		double [] setABuckets = new double[numberOfBuckets];
		double [] setBBuckets = new double[numberOfBuckets];
		
		
		
		
		
		/* create the relative expression histogram */
		double [] relativeExpression = new double[numberOfBuckets];
		for (int i = 0; i < numberOfBuckets; i++) {
			int bucketPosition = i * bucketSize;
			double bucketDifference = 0;
			if (setABuckets[i] > setBBuckets[i]) {
				if (setBBuckets[i] == 0) {
					/* if they are both zero, it is zero difference.  Otherwise, use the max bar value */
					if (setABuckets[i] == 0) {
						bucketDifference = 0;
					} else {
						bucketDifference = maxBarValue;
					}
				} else {
					bucketDifference = setABuckets[i] / setBBuckets[i];
					bucketDifference -= 1;
					if (bucketDifference > maxBarValue) bucketDifference = maxBarValue;
				}
			} else {
				if (setABuckets[i] == 0) {
					if (setBBuckets[i] == 0) {
						bucketDifference = 0;
					} else {
						bucketDifference = -maxBarValue;
					}
				} else {
					bucketDifference = setBBuckets[i] / setABuckets[i];
					bucketDifference -= 1;
					if (bucketDifference > maxBarValue) bucketDifference = maxBarValue;
					bucketDifference *= -1;
				}
			}
			
			/* creating majority bars */
			int sign = 1;
			if (bucketDifference < 0) sign = -1;
			if (Math.abs(bucketDifference) >= 0.33) {
				bucketDifference = 1;
				bucketDifference *= sign;
			} else {
				bucketDifference = 0;
			}
			
			
			relativeExpression[i] = bucketDifference;
		}
		

		
		U.p("printing results");
		try {	
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));
			
			for (int i = 0; i < numberOfBuckets; i++) {
				int bucketPosition = i * bucketSize;
				pw.println(bucketPosition + "\t" + relativeExpression[i]);
				
//				pw.println(bucketPosition + "\t" + setABuckets[i]);
			}
			
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return outputFile;
	}
	
	
	

}
