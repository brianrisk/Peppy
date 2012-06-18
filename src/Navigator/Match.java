package Navigator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import Database.Row;
import Peppy.U;

/**
 * Copyright 2012, Brian Risk
 * Released under the Netscape Public License
 * 
 * @author Brian Risk
 *
 */
public class Match extends Row implements Comparable<Match>{
	

	
	static {
		columns.put("peptideSequence", String.class);
		columns.put("spectrumFile", File.class);
		columns.put("FilePath", File.class);
		
		columns.put("matchFile", File.class);
		columns.put("spectrumMD5", String.class);
		columns.put("spectrumID", Integer.class);
		columns.put("score", Double.class);
		columns.put("amplificationScore", Double.class);
		columns.put("dataType", String.class);
		columns.put("sequenceName", String.class);
		columns.put("proteinName", String.class);
		columns.put("start", Integer.class);
		columns.put("stop", Integer.class);
		columns.put("strand", String.class);
		columns.put("isModified", Boolean.class);
		columns.put("modIndex", Integer.class);
		columns.put("modMass", Double.class);
		columns.put("modificationCount", Integer.class);
		columns.put("modificationMassArray", String.class);
		columns.put("inORF", Boolean.class);
		columns.put("sizeOfORF", Integer.class);
		columns.put("PrecursorNeutralMass", Double.class);
		columns.put("RankCount", Integer.class);
		
	}
	
	

	protected double score;
	
	public void set(String name, Object value) {
		if (name.equals("score")) {
			score = (Double) value;
		}
		super.set(name, value);
	}
	
	
	
	public int compareTo(Match other) {
		if (score > other.getScore()) return -1;
		if (score < other.getScore()) return 1;
		return 0;
	}
	
	public double getScore() {
		return score;
	}
	
	public static ArrayList<Match> loadMatches(File file) {
		ArrayList<Match> out = new ArrayList<Match>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			
			/* read the header lines */
			br.readLine();
			br.readLine();
			br.readLine();
			
			/* find out what property each column represents */
			String line = br.readLine();
			String [] propertyNames = line.split("\t");
			
			/* read in the first line */
			line = br.readLine();
			
			while (line != null) {
				String [] chunks = line.split("\t");
				Match match = new Match();
				for (int i = 0; i < propertyNames.length; i++) {
					Class<?> propertyType = columns.get(propertyNames[i]);
					if (propertyType == null  || propertyType.equals(String.class)) {
						match.set(propertyNames[i], chunks[i]);
					} else {
						try {
							if (propertyType.equals(File.class)) {
								match.set(propertyNames[i], new File(chunks[i]));
							}
							if (propertyType.equals(Integer.class)) {
								match.set(propertyNames[i], Integer.parseInt(chunks[i]));
							}
							if (propertyType.equals(Double.class)) {
								match.set(propertyNames[i], Double.parseDouble(chunks[i]));
							}
							if (propertyType.equals(Boolean.class)) {
								match.set(propertyNames[i], Boolean.parseBoolean(chunks[i]));
							}
						}
						catch (Exception e) {
							U.p("error with property " + propertyNames[i] + " in file " + file.getAbsolutePath());
						}
					}
				}
				out.add(match);
				line = br.readLine();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return out;
	}
	

	
	public boolean equals(Match other) {
		if (other.getScore() != getScore()) return false;
		if (!other.getString("spectrumMD5").equals(getString("spectrumMD5"))) return false;
		if (!other.getString("peptideSequence").equals(getString("peptideSequence"))) return false;
		return true;
	}
	

	
	
}
