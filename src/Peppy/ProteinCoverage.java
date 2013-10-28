package Peppy;

import java.util.Hashtable;

public class ProteinCoverage extends Protein implements Comparable<ProteinCoverage> {
	
	boolean [] coverage;
	double score;
	boolean identifiable = false;
	int hits = 0;
	boolean largestGapIsAtEnd = false;
	int longestGap = 0;
	int secondLongestGap = 0;
	int coveredArea = 0;
	Hashtable<String, String> uniqueFoundPeptides = new Hashtable<String, String>();
	
	public ProteinCoverage(Protein protein) {
		super(protein.getName(), protein.getAcidString(), protein.isDecoy());
		coverage = new boolean[getAcidString().length()];
	}
	
	public void addCoverage(Navigator.Match match) {
		for (int index = match.getInt("start"); index < match.getInt("stop"); index++) {
			coverage[index] = true;
		}
		hits++;
		uniqueFoundPeptides.put(match.getString("peptideSequence"), match.getString("peptideSequence"));
	}
	
	public void calculateScore() {
		/* get largest contiguous gap */
		longestGap = 0;
		secondLongestGap = 0;
		int presentGap = 0;
		int gapStart = 0;
		boolean covered;
		coveredArea = 0;
		int largestGapIndex = 0;
		
		/* find the longest gap and where it is */
		for (int index = 0; index < coverage.length; index++) {
			covered = coverage[index];
			if (covered) {
				coveredArea++;
				if (presentGap > longestGap) {
					secondLongestGap = longestGap;
					longestGap = presentGap;
					largestGapIndex = gapStart;
					if (gapStart == 0) {
						largestGapIsAtEnd = true;
					} else {
						/* have to do this as it may have been set to "true" earlier */
						largestGapIsAtEnd = false;
					}
				}
				presentGap = 0;
			} else {
				/* if we haven't measured any gap but we're in an uncovered area, then this is the start */
				if (presentGap == 0) gapStart = index;
				
				/* increment the gap size */
				presentGap ++;
			}
		}
		
		/* one final comparison */
		if (presentGap > longestGap) {
			secondLongestGap = longestGap;
			longestGap = presentGap;
			largestGapIsAtEnd = true;
		}
		
		/* find how many Rs and Ks are in protein */
		int cleavageCount = 0;
		for (int index = 0; index < this.getLength(); index++) {
			char acid = this.getAcidString().charAt(index);
			if (acid == 'R') cleavageCount++;
			if (acid == 'K') cleavageCount++;
		}
		
		int gapCleavageCount = 0;
		for (int index = largestGapIndex; index < this.getLength(); index++) {
			if (coverage[index]) break;
			char acid = this.getAcidString().charAt(index);
			if (acid == 'R') gapCleavageCount++;
			if (acid == 'K') gapCleavageCount++;
		}
		

		
		double probability = (double) (cleavageCount - gapCleavageCount + 1) / (cleavageCount + 1);
		score = Math.pow(probability, uniqueFoundPeptides.size());
		score *= (cleavageCount - gapCleavageCount + 1);
		score = - Math.log10(score);
		
		/* when there is one peptide match and it is at the very end of the protein,
		 * then it appears that the gap and the full protein have the exact same
		 * number of cleavage acids, producing a score of infinity.  Not good.
		 * adding one to the numerator removes this possiblity.  Adding one to the
		 * denominator avoids division by zero.
		 */
		
//		score = (double) longestGap * longestGap / coverage.length;
		
//		if (secondLongestGap == 0) secondLongestGap = 1;
//		score = (double) longestGap * longestGap / (secondLongestGap * secondLongestGap);
//		score = (double) coveredArea * coveredArea / coverage.length;
	}

	public int compareTo(ProteinCoverage o) {
		if (score < o.getScore()) return 1;
		if (score > o.getScore()) return -1;
		return 0;
	}

	public double getScore() {
		return score;
	}

	public boolean[] getCoverage() {
		return coverage;
	}

	public void setIdentifiable(boolean identifiable) {
		this.identifiable = identifiable;
	}

	public boolean isIdentifiable() {
		return identifiable;
	}

	public int getHits() {
		return hits;
	}

	public boolean isLargestGapIsAtEnd() {
		return largestGapIsAtEnd;
	}

	public int getLongestGap() {
		return longestGap;
	}

	public int getSecondLongestGap() {
		return secondLongestGap;
	}

}
