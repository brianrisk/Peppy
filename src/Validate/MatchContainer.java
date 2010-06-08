package Validate;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import Peppy.SpectrumPeptideMatch;
import Utilities.U;

/**
 * This is a place which will hold a SpecatrumPeptideMatch
 * as well as some meta-data such as if the match
 * is correct or not.
 * 
 * @author Brian Risk
 *
 */
public class MatchContainer implements Comparable<MatchContainer>{
	
	private boolean isTrue;
	private SpectrumPeptideMatch match;
	
	public MatchContainer(SpectrumPeptideMatch match, boolean isTrue) {
		this.isTrue = isTrue;
		this.match = match;
	}
	
	/**
	 * This determines if it is true by going finding the spectrum file
	 * for match and then locating the correct peptide sequence file.
	 * 
	 * NOTE: this method requires that the directory structure conform
	 * to specific standards.
	 * 
	 * @param match
	 */
	public MatchContainer(SpectrumPeptideMatch match) {
		//find the file for the correct peptide
		File spectrumFile = match.getSpectrum().getFile();
		File testFolder = spectrumFile.getParentFile().getParentFile();
		File peptideFolder = new File(testFolder, "peptides");
		File peptideFile = new File(peptideFolder, spectrumFile.getName());
		
		//load in the correct peptide string
		String peptideString = "";
		try {
			BufferedReader br = new BufferedReader(new FileReader(peptideFile));
			//read the first line;
			peptideString = br.readLine();
			//close;
			br.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//testing that we've got a valid peptide file
		boolean haveAnswer = true;
		if (peptideString == null) {
			haveAnswer = false;
		}
		peptideString = peptideString.trim();
		if (peptideString.equals("")) {
			haveAnswer = false;
		}
		
		//test equality
		if (haveAnswer) {
			isTrue = match.getPeptide().getAcidSequence().equals(peptideString);
		} else {
			U.p("ERROR: there was not a valid peptide file at " + peptideFile.getName());
		}
	}
	
	public double getEValue() {return match.getEValue();}
	
	public boolean isTrue() {return isTrue;}

	public int compareTo(MatchContainer o) {
		double difference = getEValue() - o.getEValue();
		//we want to sort from least to greatest
		if (difference > 0) return -1;
		if (difference < 0) return  1;
		return 0;
	}
	

}
