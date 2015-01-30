package Experimental;

/**
 * A class that stores information about fuzzy matches between a peptide
 * sequence and a protein
 * @author Brian Risk
 *
 */
public class FuzzyAlignment {
	
	String altPeptide;

	int distance; // Levenshtein distance
	int differenceIndex = -1; // if the distance is 1, then this points to the index to change to make distance 0
	
	public FuzzyAlignment(String altPeptide, int distance, int differenceIndex) {
		super();
		this.altPeptide = altPeptide;
		this.distance = distance;
		this.differenceIndex = differenceIndex;
	}

	/**
	 * 
	 * @param peptide  A string; sending in a Peptide object would be inefficient as it would have to be converted to a string each time.
	 * @param protein
	 * @param maxDistance the maximum allowed distance.  If this is exceeded, then null is returned
	 * @return
	 */
	public static FuzzyAlignment getFuzzyAlignment(String peptide, String protein, int maxDistance) {

		
		// case: 100% match
		int foundIndex = protein.indexOf(peptide);
		if (foundIndex != -1) {
			String altPeptide = protein.substring(foundIndex, foundIndex + peptide.length());
			FuzzyAlignment out = new FuzzyAlignment(altPeptide, 0, -1);
			return out;
		}
		
		//case: one acid difference
		String left;
		String right;
		int midPoint = peptide.length() / 2;
		
		// the difference is at the first acid
		right = peptide.substring(1);
		foundIndex = protein.indexOf(right);
		if (foundIndex == 0) {
			String altPeptide = "*" + right;
			FuzzyAlignment out = new FuzzyAlignment(altPeptide,1,0);
			return out;
		}
		
		// the difference is in the first half of the peptide
		for (int index = 1; index < midPoint; index++) {
			left = peptide.substring(0, index);
//			right = (index + 1, peptide.length());
		}
		
		//the difference is in the second half of the peptide
		
		//the difference is in the final acid
		
		return null;
	}



	public int getDistance() {
		return distance;
	}

	public int getDifferenceIndex() {
		return differenceIndex;
	}

}
