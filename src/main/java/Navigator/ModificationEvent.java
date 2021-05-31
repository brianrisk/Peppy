package Navigator;


/**
 * This is a modification at a specific acid on a
 * specific peptide.
 * 
 * @author Brian Risk
 *
 */
public class ModificationEvent implements Comparable<ModificationEvent>{
	
	String peptideSequence;
	int location;
	char acid;
	int modMass;
	int count = 0;
	double score;
	Match match;
	
	
	public ModificationEvent(String peptideSequence, int location, int modMass, char acid, Match match) {
		super();
		this.peptideSequence = peptideSequence;
		this.location = location;
		this.modMass = modMass;
		this.acid = acid;
		this.match = match;
	}
	
//	public void calculateScore(double mean, double sd, int numberOfModificationsOfThisType) {
//		NormalDistribution nd = new NormalDistribution(mean, sd);
//		double score = nd.density(count);
//		score *= numberOfModificationsOfThisType;
//		score = - Math.log10(score);
//	}
	
	/**
	 * 
	 * @param the number of modtypes divided by the number of peptides containing the given acid
	 */
	public void calculateScore(double probability) {
		score = Math.pow(probability, count);
		score = - Math.log10(score);
		
	}
	
	public void incrementCount() {
		count++;
	}
	
	public String getPeptideSequence() {
		return peptideSequence;
	}
	
	public int getLocation() {
		return location;
	}
	
	public int getModMass() {
		return modMass;
	}
	
	public double getScore() {
		return score;
	}

	public int compareTo(ModificationEvent arg0) {
		if (score < arg0.getScore()) return 1;
		if (score > arg0.getScore()) return -1;
		return 0;
	}

	public char getAcid() {
		return acid;
	}

	public int getCount() {
		return count;
	}

	public Match getMatch() {
		return match;
	}
	

}
