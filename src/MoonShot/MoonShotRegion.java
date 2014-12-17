package MoonShot;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;

import Navigator.Match;
import Peppy.U;

public class MoonShotRegion  {

	int newMatchLeway = 50000;
	ArrayList<Match> matches = new ArrayList<Match>();
	
	int id;
	
	int start;
	int stop;
	boolean strand;
	String sequenceName;
	String geneName;
	String interest;
	
	int maxNovelDistance = 0;
	int minUniqueCount = Integer.MAX_VALUE;
	int mostSimilarRegionID = -1;
	
	Hashtable<String, Integer> peptideTallies = new Hashtable<String, Integer>();
	Hashtable<String, Integer> cellLineTallies = new Hashtable<String, Integer>();
	Hashtable<String, Integer> compartmentTallies = new Hashtable<String, Integer>();
	Hashtable<String, Integer> fractionTallies = new Hashtable<String, Integer>();
	Hashtable<String, Match> peptides = new Hashtable<String, Match>();


	public MoonShotRegion(Match match) {

			start = match.getInt("start");
			
			stop = match.getInt("stop");
			
			//formatting sequence name (removing possible ">")
			sequenceName = match.getString("sequenceName");
			
			geneName = match.getString("geneName");
			
			interest = match.getString("interest");
			
			String strandString = match.getString("strand");
			strand = false;
			if (strandString.equals("+")) strand = true;
			
			
			addMatch(match);
	}
	
	/**
	 * Adds match if it is within the leeway window
	 * 
	 * returns true of match was added.
	 * 
	 * @param match
	 * @return
	 */
	public boolean addMatch(Match match) {
		
		boolean wasAdded = false;
		
		if (!match.getString("sequenceName").equals(sequenceName)) return false;
		
		String strandString = match.getString("strand");
		String peptide = match.getString("peptideSequence");
		int start = match.getInt("start");
		int stop = match.getInt("stop");
		
		boolean newStrand = false;
		if (strandString.equals("+")) newStrand = true;
		if (newStrand != strand) return false;
		

		if (
				match.getInt("start") >= (this.start - this.newMatchLeway) &&
				match.getInt("stop") < (this.start + this.newMatchLeway)) {
			
			//check if already added
			boolean alreadyAdded = false;
			for (Match addedMatch: matches) {
				if (addedMatch.getString("spectrumMD5").equals(match.getString("spectrumMD5"))) {
					alreadyAdded = true;
					break;
				}
			}
			
			//if not already added, add it
			if (!alreadyAdded) {
				matches.add(match);
				
				//calculate max novel distance
				int novelDistance = match.getInt("novelDistance");
				if (novelDistance > maxNovelDistance) {
					maxNovelDistance = novelDistance;
				}
				
				//tallying peptides
				Integer peptideTally = peptideTallies.get(peptide);
				if (peptideTally == null) peptideTally = 0;
				peptideTally = peptideTally + 1;
				peptideTallies.put(match.getString("peptideSequence"), peptideTally);
				
				//peptide tally is baseline for min unique
				minUniqueCount = peptideTallies.size();
				
				///sample filepath: Volumes/New Volume/1_IPAS350_H969/LO2_IP0350_Media_SG58to59.mgf
				String fileName = match.getFile("FilePath").getName();
				String [] sampleProperties = fileName.split("_");
				
				//change for moon shot [1] vs. cptac [0]
				String cellLine = sampleProperties[0];
				
				//for moon shot
//				String compartment = sampleProperties[2];
//				String fraction = sampleProperties[3];
				
				Integer cellLineTally = cellLineTallies.get(cellLine);
				if (cellLineTally == null) cellLineTally = 0;
				cellLineTally = cellLineTally + 1;
				cellLineTallies.put(cellLine, cellLineTally);
				
//				Integer compartmentTally = compartmentTallies.get(compartment);
//				if (compartmentTally == null) compartmentTally = 0;
//				compartmentTally = compartmentTally + 1;
//				compartmentTallies.put(compartment, compartmentTally);
				
				
//				Integer fractionTally = fractionTallies.get(fraction);
//				if (fractionTally == null) fractionTally = 0;
//				fractionTally = fractionTally + 1;
//				fractionTallies.put(fraction, fractionTally);
				
				if (start < this.start) this.start = start;
				if (stop > this.stop) this.stop = stop;
				
			}
			
			//regardless if we add it here, if already added still flag wasAdded
			wasAdded = true;
		}
		
		if (wasAdded) {
			peptides.put(peptide, match);
		}
		
		return wasAdded;
	}
	
	
	public String toString() {

		StringBuffer sb = new StringBuffer();
		
		//ID
		sb.append(id + "\t");
		
		//gene name
		sb.append(geneName + "\t");
		
		//interest
		//sb.append(interest + "\t");
		
		//locus: chr + start
		sb.append(sequenceName + ":" + start + "\t");

		//calculate total PSMs
		sb.append(matches.size() + "\t");
		
		//calculate total unique sequences
		sb.append(peptideTallies.size() + "\t");
		
		//maximum novel distance
		//sb.append(maxNovelDistance + "\t");
		
		//print cell line tally
		//sb.append(cellLineTallies.size() + "\t");
		
		//print min unique count
		sb.append(minUniqueCount + "\t");
		sb.append(mostSimilarRegionID + "\t");
		
		for (String peptide: peptides.keySet()) {
			sb.append(peptide + "\t");
		}
		
		return sb.toString();
	}
	
	public String toHTML() {
		StringBuffer sb = new StringBuffer();
		
		//start the row
		sb.append("<tr>");
		
		//ID
		sb.append("<td>" + id + "</td>");
		
		//gene name
		if (geneName == null) {
			sb.append("<td>" + geneName +  "</td>");
		} else {
			if (geneName.indexOf('.') != -1 || geneName.indexOf('-') != -1 || geneName.equals("null")) {
				sb.append("<td>" + geneName +  "</td>");
			} else {
				sb.append("<td><a href=\"http://www.genecards.org/cgi-bin/carddisp.pl?gene=" + geneName + "\">" + geneName + "</a></td>");
			}
		}
		
		
		//interest
		//sb.append("<td>" + interest + "</td>");
		
		//locus: chr + start
		//http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr10%3A108228413-108228413
		String formattedSequence = sequenceName;
		if (sequenceName.startsWith(">")) formattedSequence = formattedSequence.substring(1);
		String complement = "1";
		if (strand) complement = "0";
		String ucsc = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&complement_hg19=" + complement + "&position=" + formattedSequence + "%3A" + start + "-" + stop;
		sb.append("<td><a href='" + ucsc + "'>" + formattedSequence + ":" + start + "</a></td>");

		//calculate total PSMs
		sb.append("<td>" + matches.size() + "</td>");
		
		//calculate total unique sequences
		sb.append("<td>" + peptideTallies.size() + "</td>");
		
		//maximum novel distance
		//sb.append("<td>" + maxNovelDistance + "</td>");
		
		//print cell line tally
		//sb.append("<td>" + cellLineTallies.size() + "</td>");
		
		//print min unique count
		sb.append("<td>" + minUniqueCount + "</td>");
		sb.append("<td>" + mostSimilarRegionID + "</td>");
		
		//end the row
		sb.append("<tr>");
		
		return sb.toString();
	}
	
	
	public int getUniqueCount(MoonShotRegion region) {
		int nonOverlap = 0;
		for (String peptide: peptideTallies.keySet()) {
			Integer count = region.getCountForPeptide(peptide);
			if (count == null) nonOverlap++;
		}
		if (nonOverlap < minUniqueCount) {
			minUniqueCount = nonOverlap;
			mostSimilarRegionID = region.getId();
		}
		return nonOverlap;
	}
	
	public Integer getCountForPeptide(String peptide) {
		return peptideTallies.get(peptide);
	}

	public int getId() {
		return id;
	}

	public void setId(int id) {
		this.id = id;
	}
	
	public void saveHTML(File parentDirectory) {
		 try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(parentDirectory, id + ".html"))));
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public Hashtable<String, Match> getPeptides() {
		return peptides;
	}

	

}
