package Peppy;

import java.util.ArrayList;

import Utilities.U;

/**
 * Important note:  all begin indicies are inclusive, all ends are exclusive
 * @author Brian Risk
 *
 */

public class Protein implements Comparable<Protein>{
	
	private String name;
	private int start = -1;
	private String acidString;
	private boolean isSpliced;
	private int intronStart = -1;
	private int intronStop = -1;
	private boolean isForward = true;
	private ArrayList<Peptide> peptides;
	private int hitCount = 0;
	private Sequence sequence;
	
	static final int maxCleavages = Properties.numberOfMissedCleavages + 1;

	public Protein(String name, String acidString) {
		this(name, 0, acidString);
	}
	
	public Protein(String name, int start, String acidString) {
		this(name, start, acidString, false, -1, -1, true, null);
	}
	
	/**
	 * 
	 * @param name
	 * @param start
	 * @param acidString  The acid string is assumed to, at most, have only one '.' and if one does exist it can only exist at the very end of thes string.
	 * @param isSpliced
	 * @param intronStart
	 * @param intronStop
	 * @param isForward
	 * @param sequence
	 */
	public Protein(String name, int start, String acidString, boolean isSpliced, int intronStart, int intronStop, boolean isForward, Sequence sequence) {
		this.name = name;
		this.start = start;
		this.acidString = acidString;
		this.isSpliced = isSpliced;
		this.intronStart = intronStart;
		this.intronStop = intronStop;
		this.isForward = isForward;
		this.sequence = sequence;
	}
	
	public ArrayList<Peptide> digest() {
		peptides = new ArrayList<Peptide>();
		if (acidString.length() < 4) return peptides;
		
		char aminoAcid = acidString.charAt(0);
		char previousAminoAcid = aminoAcid;

		//setting up the acid indicies
		int [] acidIndicies = new int[acidString.length()];
		if (isSpliced) {
			//Setting up nucleotide positions
			//doing this as divisions at intron start/stop might not be clean; this is just the
			//easiest way to keep track of such things.
			int [] nucleotideIndicies = new int[acidString.length() * 3];
			for (int i = start; i < intronStart; i++) {
				nucleotideIndicies[i - start] = i;
			}
			for (int i = intronStop; i < acidString.length() * 3; i++) {
				nucleotideIndicies[i - start] = i;
			}
			
			if (isForward) {
				for (int i = 0; i < acidString.length(); i++) {
					acidIndicies[i] = nucleotideIndicies[i * 3];
				}
			} else {
				for (int i = acidString.length() - 1; i >= 0; i--) {
					acidIndicies[i] = nucleotideIndicies[i * 3];
				}
			}
		} else {
			for (int i = 0; i < acidString.length(); i++) {
				if (isForward) {
					acidIndicies[i] = (3 * i) + start;
				} else {
					acidIndicies[i] = start - (3 * i) - 2;
				}
			}
		}
		
		
		
		//Where we store all of our forming peptides
		ArrayList<PeptideUnderConstruction> peptidesUnderConstruction = new ArrayList<PeptideUnderConstruction>();
		
		//start the first amino acid as peptide
		peptidesUnderConstruction.add(new PeptideUnderConstruction(acidIndicies[0], aminoAcid));
		
		//special cases happen at the end of the sequence.  These indicies will come into play
		int finalIndex = acidString.length() - 1;
	
		//advance along each amino acid
		//start 1 out so we have a previous amino acid
		for (int i = 1; i < finalIndex; i++) {

			//getting the present amino acid
			aminoAcid = acidString.charAt(i);
			
			//add the present amino acid to all forming peptides
			for (PeptideUnderConstruction puc: peptidesUnderConstruction) {
				puc.addAminoAcid(aminoAcid);
			}
			
			//create new forming peptides if necessary
			if ( (isStart(aminoAcid)) ||  // start a new peptide at M
				 (isStart(previousAminoAcid) && !isStart(aminoAcid)) || // handle possible N-terminal methionine truncation products
				 (isBreak(previousAminoAcid) && !isStart(aminoAcid))  )  // Create new peptides after a break, but only if we wouldn't have created a new one with M already
			{		
				peptidesUnderConstruction.add(new PeptideUnderConstruction(acidIndicies[i], aminoAcid));
			}
			
			//if we are at a break, 
			if (isBreak(aminoAcid)) {
				if (isBreak(aminoAcid)) {
					//remove those which have exceeded the max break count
					int size = peptidesUnderConstruction.size();
					for (int pucIndex = 0; pucIndex < size; pucIndex++) {
						PeptideUnderConstruction puc = peptidesUnderConstruction.get(pucIndex);
						if (puc.getBreakCount() > maxCleavages) {
							peptidesUnderConstruction.remove(pucIndex);
							pucIndex--;
							size--;
						}
					}
				}
				
				//add forming peptides to peptide list
				for (PeptideUnderConstruction puc: peptidesUnderConstruction) {
					evaluateNewPeptide(
							puc,
							intronStart,
							intronStop,
							acidIndicies[i],
							isForward,
							name,
							peptides);
				}

			}
			

			//skip X sequences
			if (aminoAcid == 'X') {
				while (acidString.charAt(i) == 'X') {
					i++;
					if (i ==  acidString.length()) break;
				}
				i++;
			}
			
			previousAminoAcid = aminoAcid;
		}
		
		//add the final amino acid and all of the still-forming peptides
		//This handles weird edge cases
		aminoAcid = acidString.charAt(finalIndex);
		boolean addPeptide;
		for (PeptideUnderConstruction puc: peptidesUnderConstruction) {
			puc.addAminoAcid(aminoAcid);
			addPeptide = false;
			if (isStop(aminoAcid)) {
				if (isBreak(previousAminoAcid)) {
					//do nothing
				} else {
					if (puc.getBreakCount() < maxCleavages ) {
						addPeptide = true;
					}
				}
			} else {
				if (isBreak(aminoAcid)) {
					if (puc.getBreakCount() <= maxCleavages ) {
						addPeptide = true;
					}
				} else {
					if (puc.getBreakCount() < maxCleavages ) {
						addPeptide = true;
					}
				}
			}
			if (addPeptide) {
				evaluateNewPeptide(
						puc,
						intronStart,
						intronStop,
						acidIndicies[finalIndex],
						isForward,
						name,
						peptides);
			}
		}
		
		//no need for that string anymore.  Get rid of it.
		acidString = null;
		
		return peptides;
	}
	
	
	/**
	 * This first creates a peptide from the peptide under construction.
	 * This is mildly complicated as the peptide has different constructors
	 * depending on if it comes form a protein database or nucleotide sequence
	 * @param puc
	 * @param intronStartIndex
	 * @param intronStopIndex
	 * @param acidIndex
	 * @param isFromSequence
	 * @param isForward
	 * @param geneSequence
	 * @param proteinName
	 * @param peptides
	 */
	private void evaluateNewPeptide(
			PeptideUnderConstruction puc,
			int intronStartIndex,
			int intronStopIndex,
			int acidIndex,
			boolean isForward,
			String proteinName,
			ArrayList<Peptide> peptides
			) {
		boolean isSpliced = false;
		int peptideIntronStartIndex;
		int peptideIntronStopIndex;
		Peptide peptide;
		
		//splice related
		isSpliced =  (puc.getStartIndex() < intronStartIndex && acidIndex > intronStopIndex);
		if (isSpliced) {
			peptideIntronStartIndex = intronStartIndex;
			peptideIntronStopIndex = intronStopIndex;
		} else {
			peptideIntronStartIndex = -1;
			peptideIntronStopIndex = -1;
		}
		//If this is coming from DNA or RNA, there is a different peptide constructor
		peptide = new Peptide(
				puc.getSequence(),
				puc.getStartIndex(),
				acidIndex,
				peptideIntronStartIndex,
				peptideIntronStopIndex,
				isForward,
				sequence,
				isSpliced);
		
		//add peptide if it meets certain criteria
		if (peptide.getMass() >= Properties.peptideMassThreshold) {
			peptides.add(peptide);
		}
	}
	
	/**
	 * @return the peptides
	 */
	public ArrayList<Peptide> getPeptides() {
		if (peptides == null) digest();
		return peptides;
	}

	/**
	 * @return the hitCount
	 */
	public int getHitCount() {
		return hitCount;
	}
	
	public int getStart() {
		return start;
	}

	public String getAcidString() {
		return acidString;
	}

	public boolean isForward() {
		return isForward;
	}

	public int calculateHitCount() {
		hitCount = 0;
		for (Peptide peptide: peptides) {
//			if (peptide.getHitCount() > 0) hitCount++;
		}
		return hitCount;
	}

	private boolean isStart(char aminoAcid) {
		return (aminoAcid == 'M');
	}
	
	private boolean isStop(char aminoAcid) {
		return (aminoAcid == '.');
	}
	
	
	private boolean isBreak(char aminoAcid) {
		return ( aminoAcid == 'K' || aminoAcid == 'R' || aminoAcid == 'X');
	}

	public int compareTo(Protein other) {
		return other.getHitCount() - hitCount;
	}
	

}
