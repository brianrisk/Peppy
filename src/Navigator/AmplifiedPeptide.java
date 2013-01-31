package Navigator;

import Peppy.U;

/**
 * A support class for AmplificationPlot
 * 
 * @author Brian Risk
 *
 */
public class AmplifiedPeptide {
	
	String peptide;
	int count = 1;
	double score = 0;
	boolean valid = true;
	
//	public static void main(String [] args) {
//		String test = "PTBP1|PTBP1-201|531|AmpScore: 0.0|#ofVariants: 0|HasPrematureStop: false|PrematureStopLength: 0|ENST0000";
//		U.p(test.replace('|', '\t'));
//	}
	
	public AmplifiedPeptide(Match match) {
		String peptideString = match.getString("peptideSequence");
		peptide = peptideString;
		
		/* extract the amplification score from the sequence name */
		//PTBP1|PTBP1-201|531|AmpScore: 0.0|#ofVariants: 0|HasPrematureStop: false|PrematureStopLength: 0|ENST0000
		String sequenceName = match.getString("sequenceName");
		sequenceName = sequenceName.replace('|', '\t');
		String [] chunks = sequenceName.split("\t");
		if (chunks.length < 4) {
			valid = false;
		} else {
			String scoreString = chunks[3];
			if (!scoreString.startsWith("AmpScore")) {
				valid = false;
			} else {
				scoreString = scoreString.substring("AmpScore: ".length());
				score = Double.parseDouble(scoreString);
			}
			
		}
//		try {
//			
//		} catch (Exception e) {
//			U.p(chunks.length);
//			U.p(sequenceName);
//		}
		
	}
	
	
	public void increment() {
		count++;
	}
	
	public int getCount() {
		return count;
	}
	
	public double getScore() {
		return score;
	}
	
	public String getPeptide() {
		return peptide;
	}
	
	public boolean isValid() {
		return valid;
	}

}
