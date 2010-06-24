package Validate;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import Peppy.Peptide;
import Peppy.Match;
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
	private Match match;
	private String correctAcidSequence;
	private Match trueMatch = null;
	

	/**
	 * This determines if it is true by going finding the spectrum file
	 * for match and then locating the correct peptide sequence file.
	 * 
	 * NOTE: this method requires that the directory structure conform
	 * to specific standards.
	 * 
	 * @param match
	 */
	public MatchContainer(Match match) {
		this.match = match;
		//find the file for the correct peptide
		File spectrumFile = match.getSpectrum().getFile();
		File testFolder = spectrumFile.getParentFile().getParentFile();
		File peptideFolder = new File(testFolder, "peptides");
		File peptideFile = new File(peptideFolder, spectrumFile.getName());
		
		//load in the correct peptide string
		boolean validPeptideFile = true;
		String peptideString = "";
		try {
			BufferedReader br = new BufferedReader(new FileReader(peptideFile));
			//read the first line;
			peptideString = br.readLine();
			//close;
			br.close();
		} catch (FileNotFoundException e) {
			validPeptideFile = false;
			e.printStackTrace();
		} catch (IOException e) {
			validPeptideFile = false;
			e.printStackTrace();
		}
		
		//testing that we've got a valid peptide file
		if (peptideString == null) {
			validPeptideFile = false;
		}
		peptideString = peptideString.trim();
		if (peptideString.equals("")) {
			validPeptideFile = false;
		}
		
		//test equality
		if (validPeptideFile) {
			correctAcidSequence = peptideString;
			isTrue = match.getPeptide().equals(peptideString);
			
			//If this match is not the true match:
			//this is in case the peptide is not present in the database
			//also in the off case the score is equal
			if(!isTrue) {
				trueMatch = new Match(match.getSpectrum(), new Peptide(correctAcidSequence), null);
				if (trueMatch.getScoreTandemFit() == match.getScoreTandemFit()) {
					isTrue = true;
				}
			}
		} else {
			U.p("ERROR: there was not a valid peptide file at " + peptideFile.getName());
		}
	}
	
	public double getEValue() {return match.getEValue();}
	
	public Match getMatch() {return match;}
	
	public String getCorrectAcidSequence() {return correctAcidSequence;}
	
	public boolean isTrue() {return isTrue;}
	
	public Match getTrueMatch() {return trueMatch;}

	public int compareTo(MatchContainer o) {
		double difference = getEValue() - o.getEValue();
		//we want to sort from least to greatest
		if (difference > 0) return  1;
		if (difference < 0) return -1;
		return 0;
	}
	

}
