package Navigator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import Peppy.U;

public class AmplificaitonHistogram {
	
	private static int bucketSize = 2000000;
	private static double maxBarValue = 5;
	
	private static ArrayList<Match> loadMatches(String reportLocation, double scoreCutoff) {
		ArrayList<Match> matches = Match.loadMatches(new File(reportLocation));
		matches = filter(matches, scoreCutoff);
		return matches;
	}
	
	private static ArrayList<Match> filter(ArrayList<Match> matches, double scoreCutoff) {
		ArrayList<Match> out = new ArrayList<Match>(matches.size());
		for (Match match: matches) {
//			int matchCount =  match.getInt("RankCount");
//			if (matchCount != 1) continue;
			if (match.getScore() < scoreCutoff) continue;
			out.add(match);
		}
		return out;
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
	
	public static void main(String args[]) {
		doIt();
		U.p("done");
	}
	
	private static void doIt() {
		
		/* set up histogram */
		String sequenceName = "chr8";
		int chrLength = 146364022;
		int numberOfBuckets = chrLength / bucketSize + 1;
		double [] setABuckets = new double[numberOfBuckets];
		double [] setBBuckets = new double[numberOfBuckets];
		
		
		File outputFile = new File("43--WHIM16-expression.txt");
		
		U.p("loading matches");
		ArrayList<Match> matches;
		
		/* Wash U */
		matches = loadMatches("/Users/risk2/PeppyData/WashU/reports/43/WHIM16/WashU WHIM16 xeno/report.txt", 22.556025416890137);
		addMatchesToHistogram(matches, setABuckets, sequenceName);
//		matches = loadMatches("/Users/risk2/PeppyData/WashU/reports/43/WHIM2/WashU WHIM2 xeno/report.txt", 22.556025416890137);
//		addMatchesToHistogram(matches, setBBuckets, sequenceName);
		
//		matches = loadMatches("/Users/risk2/PeppyData/WashU/reports/41/WHIM16/WashU WHIM16 xeno/report.txt", 22.556025416890137);
//		addMatchesToHistogram(matches, setABuckets, sequenceName);
//		matches = loadMatches("/Users/risk2/PeppyData/WashU/reports/41/WHIM2/WashU WHIM2 xeno/report.txt", 22.556025416890137);
//		addMatchesToHistogram(matches, setBBuckets, sequenceName);
		
//		matches = loadMatches("/Users/risk2/PeppyData/WashU/reports/33/WHIM16/WashU WHIM16 xeno/report.txt", 22.556025416890137);
//		addMatchesToHistogram(matches, setABuckets, sequenceName);
//		matches = loadMatches("/Users/risk2/PeppyData/WashU/reports/33/WHIM2/WashU WHIM2 xeno/report.txt", 22.556025416890137);
//		addMatchesToHistogram(matches, setBBuckets, sequenceName);
		
		/* Johns Hopkins */
//		matches = loadMatches("/Users/risk2/PeppyData/johns-hopkins/reports/WHIM16/JHU WHIM16 hg19/report.txt", 22.98);
//		addMatchesToHistogram(matches, setABuckets, sequenceName);
//		matches = loadMatches("/Users/risk2/PeppyData/johns-hopkins/reports/WHIM2/JHU WHIM2 hg19/report.txt", 22.98);
//		addMatchesToHistogram(matches, setBBuckets, sequenceName);
		
		/* UNC */
//		matches = loadMatches("/Users/risk2/PeppyData/UNC/reports/WHIM16/UNC WHIM16 hg19/report.txt", 21.23);
//		addMatchesToHistogram(matches, setABuckets, sequenceName);
//		matches = loadMatches("/Users/risk2/PeppyData/UNC/reports/WHIM2/UNC WHIM2 hg19/report.txt", 21.23);
//		addMatchesToHistogram(matches, setBBuckets, sequenceName);
		
		/* Vanderbilt */
//		matches = loadMatches("/Users/risk2/PeppyData/vanderbilt/reports/WHIM16/Vanderbilt WHIM16 hg19/report.txt", 13.24);
//		addMatchesToHistogram(matches, setABuckets, sequenceName);
//		matches = loadMatches("/Users/risk2/PeppyData/vanderbilt/reports/WHIM2/Vanderbilt WHIM2 hg19/report.txt", 13.24);
//		addMatchesToHistogram(matches, setBBuckets, sequenceName);
		
		
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
	}
	
	
	

}
