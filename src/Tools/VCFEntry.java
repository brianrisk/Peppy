package Tools;

public class VCFEntry {
	//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	PRIMARY
	String chrom;
	int pos;
	String ref;
	String alt;
	boolean valid = true;
	
	public boolean isValid () {
		return valid;
	}
	
	public VCFEntry (String line) {
		String[] chunks = line.split("\t");
		chrom = chunks[0].trim();
		try {
			pos = Integer.parseInt(chunks[1]);
			pos++;
		} catch (NumberFormatException nfe) {
			valid = false;
		}
		ref = chunks[3].trim();

		String[] altChunks = chunks[4].split(",");
		alt = altChunks[0].trim();
		if (ref.length() != alt.length()) {valid = false;}
	}

	public String getChrom() {
		return chrom;
	}

	public int getPos() {
		return pos;
	}

	public String getRef() {
		return ref;
	}

	public String getAlt() {
		return alt;
	}
	
	public String toString() {
		return chrom + " " + pos + " " + ref + " " + alt;
	}

}
