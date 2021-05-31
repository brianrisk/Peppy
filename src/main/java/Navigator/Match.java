package Navigator;

import Database.Row;
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
		columns.put("modLocCertain", Boolean.class);
		columns.put("inORF", Boolean.class);
		columns.put("sizeOfORF", Integer.class);
		columns.put("PrecursorNeutralMass", Double.class);
		columns.put("RankCount", Integer.class);
		columns.put("novelDistance", Integer.class);
		
		
		
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
							U.p("line: " + line);
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
	
	public String toBED() {
		//chr21; strand:fwd; frame:1; start:33790123; stop:33790225
		String chr = getString("sequenceName");
		if (chr.startsWith(">")) chr = chr.substring(1);
		int start = getInt("start");
		int stop = getInt("stop");
		return chr + "\t" + start + "\t" + stop + "\t" + getString("sequence") + "\t" + Math.round(getDouble("score")) + "\t" + getString("strand");
	}
	

	
	public boolean equals(Match other) {
		if (other.getScore() != getScore()) return false;
		if (!other.getString("spectrumMD5").equals(getString("spectrumMD5"))) return false;
		if (!other.getString("peptideSequence").equals(getString("peptideSequence"))) return false;
		return true;
	}
	

	
	
}
