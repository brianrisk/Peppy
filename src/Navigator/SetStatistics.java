package Navigator;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;

import Peppy.Peptide;
import Peppy.U;


/**
 * A quick little class to answer some questions for Karen Anderson and David Fenyo
 * 
 * Copyright 2012, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class SetStatistics {
	
	public static void main(String arg[]) {
//		massErrorReport();
		
		U.p("done");
	}
	
	
	
	public static void massErrorReport() {
		U.p("loading results");

//		double zeroFDR = 24.076968727972464; // mayo
		double zeroFDR = 20.076968727972464; // yale

//		ArrayList<Match> hg19 = Match.loadMatches(new File("/Users/risk2/PeppyData/WashU/reports/WHIM16/WashU WHIM16 hg19/report.txt"));
//		ArrayList<Match> hg19 = Match.loadMatches(new File("/Users/risk2/PeppyData/Mayo/reports/split_1331843955197/report.txt"));
		ArrayList<Match> hg19 = Match.loadMatches(new File("/Users/risk2/PeppyData/yale/reports/yale-genome 2011-11/yale-genome 2011-11_report.txt"));
		
		
		ArrayList<Double> differences = new ArrayList<Double>();
		
		for (Match match: hg19) {
			if (match.getScore() > zeroFDR && !match.getBoolean("isModified")) {
				Peptide peptide = new Peptide(match.getString("peptideSequence"));
	
				double precursor = match.getDouble("PrecursorNeutralMass");
				double difference = peptide.getMass() - precursor;
				double ppm = (difference / precursor) * 1000000;
				differences.add(ppm);
			}

		}
		Collections.sort(differences);
		
		double average = 0;
		for (double difference: differences) {
			average += difference;
		}
		U.p("min: " + differences.get(0));
		U.p("max:" + differences.get(differences.size() - 1));
		average /= differences.size();
		
		
		double sd = 0;
		for (double difference: differences) {
			double delta = average - difference;
			sd += delta * delta;
		}
		sd /= differences.size();
		sd = Math.sqrt(sd);
		
		U.p("total at this threshold: " + differences.size());
		U.p("median: " + differences.get(differences.size() / 2));
		U.p("mean: " + average);
		U.p("standard deviation: " + sd);
		
//		try {
//			PrintWriter massErrorsList = new PrintWriter(new BufferedWriter(new FileWriter("mass errors list.txt")));
//			for (double difference: differences) {
//				massErrorsList.println(difference);
//			}
//			massErrorsList.flush();
//			massErrorsList.close();
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
		
	}
	
	public static void orfStats() {
		U.p("loading results");
		

//		ArrayList<Match> hg19 = Match.loadMatches(new File("/Users/risk2/PeppyData/WashU/reports/WHIM16/WashU WHIM16 hg19/report.txt"));
		ArrayList<Match> hg19 = Match.loadMatches(new File("/Users/risk2/PeppyData/yale/reports/yale-genome 2011-11/yale-genome 2011-11_report.txt"));
		
		int inORFCount = 0;
		double ORFSize = 0;
		
		ArrayList<Integer> orfs = new ArrayList<Integer>();
		
		for (Match match: hg19) {
			if (match.getBoolean("inORF")) {
				inORFCount++;
				int size = match.getInt("SizeOfORF");
				ORFSize += size;
				orfs.add(size);
			}
		}
		Collections.sort(orfs);
		
		double averageORF = ORFSize / inORFCount;
		
		double sd = 0;
		for (Match match: hg19) if (match.getBoolean("inORF")) {
			int size = match.getInt("SizeOfORF");
			double difference = averageORF - size;
			sd += difference * difference;
		}
		sd /= inORFCount;
		sd = Math.sqrt(sd);
		
		U.p("median ORF: " + orfs.get(orfs.size() / 2));
		U.p("average ORF: " + averageORF);
		U.p("standard deviation: " + sd);
		U.p("proportion in an ORF: " + ((double) inORFCount / hg19.size()));
		U.p("total in ORF: " + (inORFCount));
		U.p("total matches: " + (hg19.size()));
	}

}
