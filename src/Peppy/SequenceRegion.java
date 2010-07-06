package Peppy;

public class SequenceRegion implements Comparable<SequenceRegion>{
	
	int startIndex;
	int numberOfNucleotides;
	int numberOfHits;
	double score;
	int [] hitHistogram;
	
	public SequenceRegion(int startIndex, int numberOfNucleotides) {
		this.startIndex = startIndex;
		this.numberOfNucleotides = numberOfNucleotides;
		hitHistogram = new int[numberOfNucleotides];
		for (int i = 0; i < numberOfNucleotides; i++) {
			hitHistogram[i] = 0;
		}
	}
	
	public void addHit(Match match) {
		int matchIndex = match.getPeptide().getIndex();
		
		//multiply by 3 because there are three nucleotides for each amino acid
		int peptideLength = match.getPeptide().getAcidSequence().length() * 3;
		
		double hitValue = Math.abs(Math.log(match.getEValue()));
		
		score += hitValue;
		
		int matchStart =  matchIndex - startIndex;
		int matchEnd = matchStart + peptideLength;
		for (int i = matchStart; i < matchEnd; i++) {
			hitHistogram[i]++;
		}
	}
	
	public int [] getHistogram() {return hitHistogram;}
	public double getScore() {return score;}

	public int compareTo(SequenceRegion o) {
		//want to sort from greatest to least
		if (score > o.getScore()) return -1;
		if (score < o.getScore()) return  1;
		return 0;
	}
	

}
