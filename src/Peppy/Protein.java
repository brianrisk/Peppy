package Peppy;

import java.util.ArrayList;

/**
 * A protein is in charge of digestion and holding its digested peptides.
 * 
 * Important note:  all begin indicies are inclusive, all ends are exclusive
 * 
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */

public class Protein {
	
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
	private boolean isDecoy;
	private double averageMass = -1;

	
	
	static final int maxCleavages = Properties.numberOfMissedCleavages + 1;

	public Protein(String name, String acidString, boolean isDecoy) {
		this(name, 0, acidString, false, -1, -1, true, null, isDecoy);
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
	public Protein(String name, int start, String acidString, boolean isSpliced, int intronStart, int intronStop, boolean isForward, Sequence_DNA sequence_DNA, boolean isDecoy) {
		this.name = name;
		this.start = start;
		this.acidString = acidString;
		this.isSpliced = isSpliced;
		this.intronStart = intronStart;
		this.intronStop = intronStop;
		this.isForward = isForward;
		this.sequence_DNA = sequence_DNA;
		this.isDecoy = isDecoy;
	}
	
	
	private ArrayList<Peptide> digest() {
		peptides = new ArrayList<Peptide>();
		
		//if our peptide has only 3 acids, return an empty list
		if (acidString.length() < 4) return peptides;
		
		/* track if in ORF */
		boolean inORF = false;
		int ORFSize = 0;
		
		char aminoAcid = acidString.charAt(0);
		char nextAminoAcid = acidString.charAt(1);
		char previousAminoAcid = '.';
		if (aminoAcid == 'M') {
			inORF = true;
			ORFSize = acidString.length();
		}

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
		peptidesUnderConstruction.add(new PeptideUnderConstruction(acidIndicies[0], aminoAcid, nextAminoAcid, inORF, ORFSize, previousAminoAcid));
		
		/* special cases happen at the end of the sequence. 
		 * Don't worry, the final amino acid will eventually be added!
		 */
		int finalIndex = acidString.length() - 1;
		
		/* set our previous amino acid */
		previousAminoAcid = acidString.charAt(0);
	
		//advance along each amino acid
		//start 1 out so we have a previous amino acid
		for (int i = 1; i < finalIndex; i++) {

			//getting the present amino acid
			aminoAcid = acidString.charAt(i);
			nextAminoAcid = acidString.charAt(i + 1);
			
			
			/* see if we are in ORF */
			if (aminoAcid == 'M') {
				inORF = true;
				ORFSize = acidString.length() - i;
			}
			
			//add the present amino acid to all forming peptides
			for (PeptideUnderConstruction puc: peptidesUnderConstruction) {
				puc.addAminoAcid(aminoAcid, nextAminoAcid);
			}
			
			//create new forming peptides if necessary
//			if (Properties.isSequenceFileDNA) {
//				if ( (isStart(aminoAcid)) ||  // start a new peptide at M
//					 (isStart(previousAminoAcid) && !isStart(aminoAcid)) || // handle possible N-terminal methionine truncation products
//					 (isBreak(previousAminoAcid, nextAminoAcid) && !isStart(aminoAcid))  )  // Create new peptides after a break, but only if we wouldn't have created a new one with M already
//				{		
//					peptidesUnderConstruction.add(new PeptideUnderConstruction(acidIndicies[i], aminoAcid, nextAminoAcid, inORF, ORFSize, previousAminoAcid));
//				}
//
//				
//			/* 'M' Does not mean a new peptide should form in proteins 
//			 * Also, we don't need to explicitly say start with 'M' because in the lines above we have that the beginning of every protein begins a new peptide */
//			} else {
				// Create new peptides after a break
				if (isBreak(previousAminoAcid, nextAminoAcid)) {		
					peptidesUnderConstruction.add(new PeptideUnderConstruction(acidIndicies[i], aminoAcid, nextAminoAcid, inORF, ORFSize, previousAminoAcid));
				} else {
					/* this is for when, in proteins, sometimes the first M is dropped */
					if (i == 1 && isStart(previousAminoAcid)) {
						peptidesUnderConstruction.add(new PeptideUnderConstruction(acidIndicies[i], aminoAcid, nextAminoAcid, inORF, ORFSize, previousAminoAcid));
					}
				}
//			}
			
			//if we are at a break, 
			if (isBreak(aminoAcid, nextAminoAcid)) {
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
			puc.addAminoAcid(aminoAcid, '.');
			addPeptide = false;
			if (isStop(aminoAcid)) {
				if (isBreak(previousAminoAcid, aminoAcid)) {
					//do nothing
				} else {
					if (puc.getBreakCount() < maxCleavages ) {
						addPeptide = true;
					}
				}
			} else {
				if (isBreak(aminoAcid, '.')) {
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
		
		/* mark all peptides as coming from a decoy database, if that is the case */
		if (isDecoy) {
			for (Peptide peptide: peptides) {
				peptide.setDecoy(true);
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
		
		/*
		 * If we are searching for selenocysteine, we are only
		 * allowing peptides with those peptide for the
		 * DNA searches
		 */
		if (Properties.useSelenocysteine && Properties.isSequenceFileDNA) {
			if (sequenceString.indexOf('U') == -1) return;
		}
		
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
			
			int stopIndex = acidIndex;
			int intronStop = peptideIntronStopIndex;
			
			/* Adjust indicies for DNA and RNA */
			if (Properties.isSequenceFileDNA) {
				if (isForward) {
					//2, not 3, because index already at 1
					stopIndex += 2;
					intronStop += 2;
				} else {
					//4 because we're already at the 1 index of the next peptide
					stopIndex -= 4;
					intronStop -= 4;
				}
			}
			
			/* create the peptide */
			peptide = new Peptide(
					sequenceString,
					puc.getStartIndex(),
					stopIndex,
					peptideIntronStartIndex,
					intronStop,
					isForward,
					sequence_DNA,
					this,
					isSpliced,
					puc.isInORF(),
					puc.getORFSize(),
					puc.getPreviousAminoAcid()
					);
			
			//add peptide if it meets certain criteria
			if (peptide.getMass() >= Properties.peptideMassMinimum && peptide.getMass() <= Properties.peptideMassMaximum) {
				peptides.add(peptide);
			}
		}
	}
	
	public String getAcidString() {
		if (acidString != null) return acidString;
		return AminoAcids.getStringForByteArray(acidByteArray);
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


	public int getStart() {
		return start;
	}
	
	public static boolean isBreak(char aminoAcid, char nextAminoAcid) {
		if (Properties.cleavageAtCarboxylSide) {
			for (char cleavageAcid: Properties.cleavageAcidList) {
				if (aminoAcid == cleavageAcid) return true;
			}
		}
		return false;
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
	
	
	
	public double getAverageMass() {
		if (averageMass <0 ) calculateAverageMass();
		return averageMass;
		
	}
	public double calculateAverageMass() {
		averageMass = 0;
		for (int index = 0; index < acidString.length(); index++) {
			averageMass += AminoAcids.getWeightAverage(acidString.charAt(index));
		}
		return averageMass;
	}


	public boolean isDecoy() {
		return isDecoy;
	}
	
	public int getLength() {
		return acidString.length();
	}
	

}
