package Peppy;

import java.util.ArrayList;

/**
 * Important note:  all begin indicies are inclusive, all ends are exclusive
 * @author Brian Risk
 *
 */

public class Protein implements Comparable<Protein>{
	
	String name;
	int start = -1;
	String acidString;
	boolean isSpliced;
	int intronStart = -1;
	int intronStop = -1;
	int intronLength = -1;
	boolean isForward = true;
	ArrayList<Peptide> peptides;
	int hitCount = 0;
	Sequence sequence;
	
	public Protein(String name, int start, String acidString) {
		this(name, start, acidString, false, -1, -1, true, null);
	}
	
	public Protein(String name, int start, String acidString, boolean isSpliced, int intronStart, int intronStop, boolean isForward, Sequence sequence) {
		this.name = name;
		this.start = start;
		this.acidString = acidString;
		this.isSpliced = isSpliced;
		this.intronStart = intronStart;
		this.intronStop = intronStop;
		intronLength = intronStop - intronStart;
		this.isForward = isForward;
		this.sequence = sequence;
	}
	
	public ArrayList<Peptide> digest() {
		ArrayList<Peptide> peptides = new ArrayList<Peptide>();
		
		char aminoAcid = acidString.charAt(0);
		char previousAminoAcid = aminoAcid;

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
				acidIndicies[i] = (i * 3) + start;
			}
		}
		
		//setting up the acid indicies
		
		
		
		
		//Where we store all of our forming peptides
		ArrayList<PeptideUnderConstruction> peptidesUnderConstruction = new ArrayList<PeptideUnderConstruction>();
		
		//start the first amino acid as peptide
		peptidesUnderConstruction.add(new PeptideUnderConstruction(acidIndicies[0], aminoAcid));
	
		//advance along each amino acid
		//start 1 out so we have a previous amino acid
		for (int i = 1; i < acidString.length(); i++) {

			//getting the present amino acid
			aminoAcid = acidString.charAt(i);
			
			//add the present amino acid to all forming peptides
			for (PeptideUnderConstruction puc: peptidesUnderConstruction) {
				puc.addAminoAcid(aminoAcid);
			}
			
			//create new forming peptides if necessary
			if ((isStart(aminoAcid)) ||  // start a new peptide at M
				(isStart(previousAminoAcid) && !isStart(aminoAcid)) || // handle possible N-terminal methionine truncation products
				(isBreak(previousAminoAcid) && !isStart(aminoAcid)))  // Create new peptides after a break, but only if we wouldn't have created a new one with M already
			{		
				peptidesUnderConstruction.add(new PeptideUnderConstruction(acidIndicies[i], aminoAcid));
			}
			

			
			//if we are at a break, add forming peptides to peptide list
			else {
				if (isBreak(aminoAcid)) {
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
			}
			
			//if stop, then clear list of peptides under construction
			if (isStop(aminoAcid)) {
				peptidesUnderConstruction = new ArrayList<PeptideUnderConstruction>();
			}
			
			//remove all peptide under construction that have reached their maximum break count	
			else {
				int size = peptidesUnderConstruction.size();
				for (int pucIndex = 0; pucIndex < size; pucIndex++) {
					PeptideUnderConstruction puc = peptidesUnderConstruction.get(pucIndex);
					if (puc.getBreakCount() > Properties.numberOfMissedCleavages) {
						peptidesUnderConstruction.remove(pucIndex);
						pucIndex--;
						size--;
					}
				}
			}
			
			//skip X sequences
			if (acidString.charAt(i) == 'X') {
				while (acidString.charAt(i) == 'X') {
					i++;
					if (i ==  acidString.length()) break;
				}
				i++;
			}
		}
		
		//adding all the remaining peptides under construction
		for (PeptideUnderConstruction puc: peptidesUnderConstruction) {
			evaluateNewPeptide(
					puc,
					intronStart,
					intronStop,
					acidIndicies[acidIndicies.length - 1],
					isForward,
					name,
					peptides);
		}
		
		return peptides;
	}
	
	public ArrayList<Peptide> reverseDigest() {
		StringBuffer reverseBuffer = new StringBuffer();
		for (int i = acidString.length() - 1; i >=0; i--) {
			reverseBuffer.append(acidString.charAt(i));
		}
		acidString = reverseBuffer.toString();
		return digest();
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
		//public Peptide(String acidSequence, int startIndex, int stopIndex, int intronStartIndex, int intronStopIndex, boolean forward, Sequence parentSequence, boolean isSpliced) {
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
		return peptides;
	}

	/**
	 * @return the hitCount
	 */
	public int getHitCount() {
		return hitCount;
	}
	
	public int calculateHitCount() {
		hitCount = 0;
		for (Peptide peptide: peptides) {
//			if (peptide.getHitCount() > 0) hitCount++;
		}
		return hitCount;
	}

	private static boolean isStart(char aminoAcid) {
		return (aminoAcid == 'M');
	}
	
	private static boolean isStop(char aminoAcid) {
		return (aminoAcid == '.');
	}
	
	private static boolean isBreak(char aminoAcid) {
		return (aminoAcid == '.' || aminoAcid == 'K' || aminoAcid == 'R' || aminoAcid == 'X');
	}

	public int compareTo(Protein other) {
		return other.getHitCount() - hitCount;
	}
	

}
