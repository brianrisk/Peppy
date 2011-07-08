package Peppy;

import java.util.ArrayList;

/**
 * A protein is in charge of digestion and holding its digested peptides.
 * 
 * Important note:  all begin indicies are inclusive, all ends are exclusive
 * 
 * @author Brian Risk
 *
 */

public class Protein implements Comparable<Protein>{
	
	private String name;
	private int start = -1;
	private String acidString;
	private byte [] acidByteArray = null;
	private boolean isSpliced;
	private int intronStart = -1;
	private int intronStop = -1;
	private boolean isForward = true;
	private ArrayList<Peptide> peptides;
	private Sequence_DNA sequence_DNA;
	private ArrayList<Match> matchesAll = new ArrayList<Match>();
	private ArrayList<Match_IMP_VariMod> matchesWithModifications = new ArrayList<Match_IMP_VariMod>();
	private ArrayList<Match> matchesWithoutModifications = new ArrayList<Match>();
	private double score = 0;
	private int [] matchPositions = null;
	private double matchCoverage = -1;
	private int matchArea = -1;
	private double modCoverage = -1;
	private int modArea = -1;
	
	private static int tracker = 0;
	public static final int T_NOTHING= tracker++;
	public static final int T_FPRXX = tracker++;
	public static final int T_FPR01 = tracker++;
	public static final int T_FPR05 = tracker++;
	public static final int T_MOD = tracker++;
	
	
	static final int maxCleavages = Properties.numberOfMissedCleavages + 1;

	public Protein(String name, String acidString) {
		this(name, 0, acidString, false, -1, -1, true, null);
	}


	/**
	 * 
	 * @param name
	 * @param start
	 * @param acidString  The acid string is assumed to, at most, have only one '.' and if one does exist it can only exist at the very end of the string.
	 * @param isSpliced
	 * @param intronStart
	 * @param intronStop
	 * @param isForward
	 * @param sequence_DNA
	 */
	public Protein(String name, int start, String acidString, boolean isSpliced, int intronStart, int intronStop, boolean isForward, Sequence_DNA sequence_DNA) {
		this.name = name;
		this.start = start;
		this.acidString = acidString;
		this.isSpliced = isSpliced;
		this.intronStart = intronStart;
		this.intronStop = intronStop;
		this.isForward = isForward;
		this.sequence_DNA = sequence_DNA;
	}
	
	/**
	 * Adds the match to our list of matches.  Updates our score;
	 * @param match
	 */
	public void addMatch(Match match) {
		matchesAll.add(match);
		matchesWithoutModifications.add(match);
		//minus because this will be negative for good e values
//		score -= Math.log(match.getEValue());

		if (matchPositions == null) {
			matchPositions = new int[acidByteArray.length];
			for (int i = 0; i < matchPositions.length; i++) {
				matchPositions[i] = T_NOTHING;
			}
		}
		int type = T_FPRXX;
		//TODO do not hard code these FPRs!
		if (match.getEValue() < 0.03359587957603186) type = T_FPR05;
		if (match.getEValue() < 0.0034847330927202927) type = T_FPR01;
		for (int i =  match.getPeptide().getStartIndex(); i < match.getPeptide().getStopIndex(); i++) {
			matchPositions[i] = type;
		}
		
		score += 1;
	}
	
	public void addMatchPTM(Match_IMP_VariMod match) {
		matchesAll.add(match);
		matchesWithModifications.add(match);
		//minus because this will be negative for good e values
//		score -= Math.log(match.getEValue());


		if (matchPositions == null) {
			matchPositions = new int[acidByteArray.length];
			for (int i = 0; i < matchPositions.length; i++) {
				matchPositions[i] = T_NOTHING;
			}
		}
		for (int i =  match.getPeptide().getStartIndex(); i < match.getPeptide().getStopIndex(); i++) {
			if (matchPositions[i] != T_FPR05 && matchPositions[i] != T_FPR01 )
			matchPositions[i] = T_MOD;
		}
		
		score += 1;
	}
	
	public ArrayList<Peptide> getUnfoundPeptides() {
		if (matchPositions != null) {
			ArrayList<Peptide> unfoundPeptides = new ArrayList<Peptide>();
			boolean peptideFound;
			for (Peptide peptide: peptides) {
				peptideFound = true;
	
				//this says that if there is any amino acid that hasn't been accounted for
				//that exists in this peptide, then the peptide hasn't been matched
				for (int i =  peptide.getStartIndex(); i < peptide.getStopIndex(); i++) {
					if (matchPositions[i] == T_NOTHING) {
						peptideFound = false;
						break;
					}
				}
				if (!peptideFound) {
					unfoundPeptides.add(peptide);
				}
			}
			return unfoundPeptides;
		}
		return peptides;
	}
	
	public double getMatchCoverage() {
		if (matchCoverage < 0) {
			matchCoverage = (double) getMatchArea() / matchPositions.length;
		}
		return matchCoverage;
	}
	
	public int getMatchArea() {
		if (matchArea < 0) {
			matchArea = 0;
			for (int i = 0; i < matchPositions.length; i++) {
				if (matchPositions[i] == T_FPR01) matchArea++;
				if (matchPositions[i] == T_FPR05) matchArea++;
			}
		}
		return matchArea;
	}
	
	public double getModCoverage() {
		if (modCoverage < 0) {
			modCoverage = (double) getModArea() / matchPositions.length;
		}
		return modCoverage;
	}
	
	public int getModArea() {
		if (modArea < 0) {
			modArea = 0;
			for (int i = 0; i < matchPositions.length; i++) {
				if (matchPositions[i] == T_MOD) modArea++;
			}
		}
		return modArea;
	}
	
	public int [] getMatchPositions() {return matchPositions;}
	
	
	public int compareTo(Protein other) {
		if (other.getScore() < getScore()) return -1;
		if (other.getScore() > getScore()) return  1;
		return 0;
	}
	
	private ArrayList<Peptide> digest() {
		peptides = new ArrayList<Peptide>();
		
		//if our peptide has only 3 acids, return an empty list
		if (acidString.length() < 4) return peptides;
		
		char aminoAcid = acidString.charAt(0);
		char previousAminoAcid = aminoAcid;

		//setting up the acid indicies
		int [] acidIndicies = new int[acidString.length()];
		if (Properties.isSequenceFileDNA) {
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
		} else {
			for (int i = 0; i < acidString.length(); i++) {
				acidIndicies[i] = i;
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
			if (Properties.isSequenceFileDNA) {
				if ( (isStart(aminoAcid)) ||  // start a new peptide at M
					 (isStart(previousAminoAcid) && !isStart(aminoAcid)) || // handle possible N-terminal methionine truncation products
					 (isBreak(previousAminoAcid) && !isStart(aminoAcid))  )  // Create new peptides after a break, but only if we wouldn't have created a new one with M already
				{		
					peptidesUnderConstruction.add(new PeptideUnderConstruction(acidIndicies[i], aminoAcid));
				}
				
			/* 'M' Does not mean a new peptide should form in proteins */
			} else {
				// Create new peptides after a break
				if (isBreak(previousAminoAcid)) {		
					peptidesUnderConstruction.add(new PeptideUnderConstruction(acidIndicies[i], aminoAcid));
				}
			}
			
			//if we are at a break, 
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
				
				//add forming peptides to peptide list
				for (PeptideUnderConstruction puc: peptidesUnderConstruction) {
					evaluateNewPeptide(
							puc,
							intronStart,
							intronStop,
							acidIndicies[i] + 1,
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
						acidIndicies[finalIndex] + 1,
						isForward,
						name,
						peptides);
			}
		}
		
		//no need for that string anymore.  Get rid of it.
		acidByteArray = AminoAcids.getByteArrayForString(acidString);
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
		String sequenceString = puc.getSequence();
		
		//splice related
		isSpliced =  (puc.getStartIndex() < intronStartIndex && acidIndex > intronStopIndex);
		if (isSpliced) {
			peptideIntronStartIndex = intronStartIndex;
			peptideIntronStopIndex = intronStopIndex;
		} else {
			peptideIntronStartIndex = -1;
			peptideIntronStopIndex = -1;
		}
		
		//might not be considering very long peptides
		if (sequenceString.length() <= Properties.maxPeptideLength) {
			//If this is coming from DNA or RNA, there is a different peptide constructor
			int stopIndex = acidIndex;
			int intronStop = peptideIntronStopIndex;
			if (isForward) {
				//2, not 3, because index already at 1
				stopIndex += 2;
				intronStop += 2;
			} else {
				//4 because we're already at the 1 index of the next peptide
				stopIndex -= 4;
				intronStop -= 4;
			}
			peptide = new Peptide(
					sequenceString,
					puc.getStartIndex(),
					stopIndex,
					peptideIntronStartIndex,
					intronStop,
					isForward,
					sequence_DNA,
					this,
					isSpliced);
			
			//add peptide if it meets certain criteria
			if (peptide.getMass() >= Properties.peptideMassThreshold) {
				peptides.add(peptide);
			}
		}
	}
	
	public String getAcidString() {
		return AminoAcids.getStringForByteArray(acidByteArray);
	}


	public double getScore() {
		if (matchPositions == null) return 0.0;
		return (double) getMatchArea() * getMatchArea() / acidByteArray.length;
//		return score;
	}

	public String getName() {
		return name;
	}

	/**
	 * @return the peptides
	 */
	public ArrayList<Peptide> getPeptides() {
		if (peptides == null) digest();
		return peptides;
	}

	public ArrayList<Match> getMatchesAll() {
		return matchesAll;
	}


	public ArrayList<Match_IMP_VariMod> getMatchesWithModifications() {
		return matchesWithModifications;
	}


	public ArrayList<Match> getMatchesWithoutModifications() {
		return matchesWithoutModifications;
	}


	public int getStart() {
		return start;
	}
	
	public static boolean isBreak(char aminoAcid) {
		return ( aminoAcid == 'K' || aminoAcid == 'R' || aminoAcid == 'X');
//		return ( aminoAcid == 'Y' || aminoAcid == 'W' || aminoAcid == 'F' || aminoAcid == 'X');
	}
	
	
	public boolean isForward() {
		return isForward;
	}

	private boolean isStart(char aminoAcid) {
		return (aminoAcid == 'M');
	}

	private boolean isStop(char aminoAcid) {
		return (aminoAcid == '.');
	}
	

}
