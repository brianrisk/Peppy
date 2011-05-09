package Validate;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import Peppy.Match;
import Peppy.Peptide;
import Peppy.Properties;
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
		String correctAcidSequence = "";
		try {
			BufferedReader br = new BufferedReader(new FileReader(peptideFile));
			//read the first line;
			correctAcidSequence = br.readLine();
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
		if (correctAcidSequence == null) {
			validPeptideFile = false;
		}
		correctAcidSequence = correctAcidSequence.trim();
		if (correctAcidSequence.equals("")) {
			validPeptideFile = false;
		}
		
		//test equality
		if (validPeptideFile) {
			this.correctAcidSequence = correctAcidSequence;
			isTrue = match.getPeptide().equals(correctAcidSequence);
			
			//If this match is not the true match:
//			if(!isTrue) {
//				trueMatch = Properties.matchConstructor.createMatch(match.getSpectrum(), new Peptide(this.correctAcidSequence));
//				double trueEValue = match.getSpectrum().getEValue(trueMatch.getScore());
//				trueMatch.setEValue(trueEValue);
//			}
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
		return match.compareTo(o.getMatch());
	}
	

}
