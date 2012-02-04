package Navigator;

public class Match_fromResults implements Comparable<Match_fromResults> {
	
	
	public static final int TYPE_PROTEIN = 0;
	public static final int TYPE_GENOME = 1;
	
	private String peptide;
	private String spectrumMD5;
	private String spectrumFile;
	private double score;
	private int matchType;
	
	
	public Match_fromResults(String peptide, String spectrumMD5, String spectrumFile,
			double score, int matchType) {
		super();
		this.peptide = peptide;
		this.spectrumMD5 = spectrumMD5;
		this.spectrumFile = spectrumFile;
		this.score = score;
		this.matchType = matchType;
	}


	public String getPeptide() {
		return peptide;
	}


	public String getSpectrumMD5() {
		return spectrumMD5;
	}


	public String getSpectrumFile() {
		return spectrumFile;
	}


	public double getScore() {
		return score;
	}


	public int getMatchType() {
		return matchType;
	}


	public int compareTo(Match_fromResults arg0) {
		if (score > arg0.getScore()) return -1;
		if (score < arg0.getScore()) return 1;
		return 0;
	}
	
	public String toString() {
		StringBuffer out = new StringBuffer();
		out.append(peptide);
		out.append('\t');
		if (matchType == TYPE_PROTEIN) {
			out.append("protein");
		} else {
			out.append("genome");
		}
		out.append('\t');
		out.append(score);
		out.append('\t');
		out.append(spectrumFile);

		
		return out.toString();
	}


}
