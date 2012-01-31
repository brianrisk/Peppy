package HMMScore;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Vector;

import Peppy.AminoAcids;
import Peppy.Peak;
import Peppy.Properties;
import Peppy.Spectrum;

public class HMMClass {
	private static int numberOfState = Defines.numberOfIons + 1;
	private static double[][] logA = new double[numberOfState][numberOfState];
	private static double[][] logEM = new double[numberOfState][Defines.bin_count];
	private static double[][] logEI = new double[numberOfState][Defines.bin_count];
	private static double[] logNTerCleavageProb = new double[200];
	private static double[] logCTerCleavageProb = new double[200];
	private static final float ZERO_FINAL = (float) 1.0e-04;
	private static float[] probArray = new float[numberOfState];
	private static float[] observedProbArray = new float[numberOfState];
	private static double[] logObservedProb = new double[numberOfState];
	
	private String acidSequence;
	private Spectrum spectrum;
	private Vector<Double> fragsMasses = new Vector<Double>();
	private ArrayList<IonMatch> fragsMatches = new ArrayList<IonMatch>();
	

	public HMMClass(String acidSequence, Spectrum spectrum) {
		//check the last character and if that is a stop, then ignore that last character
		if (acidSequence.endsWith(".") ) acidSequence = acidSequence.substring(0, acidSequence.length()-1);
//System.out.println("NOW  "+acidSequence);
		this.acidSequence = acidSequence;	
		this.spectrum = spectrum;
	}

	/***********************************************************************************************************/
	public static void calculateAve() {
		float totalNTerCleavage = 0, totalCTerCleavage = 0;
		for (int j = 0; j < Defines.numberOfAminoAcid; j++) {
			totalNTerCleavage += logNTerCleavageProb[j];
			totalCTerCleavage += logCTerCleavageProb[j];
		}
	}

	/***************************************************************************************************/
	public static int calculateInternalCount(String sequence) {

		int length = sequence.length();
		int internalIonCount = 0;
		for (int i = 1; i < length - 2; i++) {
			internalIonCount += i;
		}
		return internalIonCount;
	}

	/*****************************************************************************************************************
	 * This method does the Viterbi calculation and calculate maximum
	 * probability
	 ******************************************************************************************************************/
	
	public double score() {
	
		int i, k, kk, l, kmax;
	//	int internalMatchCount;
	//	int matchedPeakCount = 0;
		double prod, maxprod = 0.0, currentMassEmitProb, currentIntensityEmitProb;
		double currentValue, thisValue;
		int currentMassBin = 0, currentIntensityBin = 0;
		int currentIonType, currentIonIndex;
		double currentMass, currentIntensity;
		double NTerCleaveProb = 0.0, CTerCleaveProb = 0.0;
		double precursorMass = spectrum.getMass(); // (Float)
		// tandemData.get("PrecursorMass");
		ArrayList<Peak> peaks = spectrum.getPeaks();
	
//		peaks = EmitionProb.cleanPeakF	orHighMass(peaks, precursorMass);
//		if (peaks.size() > Defines.TOTAL_PEAK_COUNT) {
//			peaks = EmitionProb.getHighIntensityPeaks(peaks);
//		}
		int massSize = peaks.size();
	
		int[][] traceBack = new int[massSize + 1][numberOfState];
		double[][] viterbi = new double[numberOfState][massSize + 1];
		int[] path = new int[massSize]; // double toAdd;
		int[] ionTypeCount = new int[numberOfState - 1];
	
		for (i = 0; i <= massSize; i++) {
			for (l = 0; l < numberOfState; l++) {
				traceBack[i][l] = 0;
				viterbi[l][i] = Double.NEGATIVE_INFINITY;
			}
		}
		viterbi[0][0] = 0; // log(1) =0
		for (k = 1; k < numberOfState; k++) {
			viterbi[k][0] = Double.NEGATIVE_INFINITY; // log (0) is - inf
		}
	
		for (i = 1; i <= massSize; i++) {
			viterbi[0][i] = Double.NEGATIVE_INFINITY; // log (0) is - inf
		}
	
		generateFragmentsForSequence(acidSequence, 0);
	
		Vector<Double> sortedIntensity = EmitionProb.sortIntensities(peaks);
		 double totalIntensity = getTotalIntensity (peaks);
		int[] prevPath = { 0, 0 }, currentPath = { 0, 0 };
		int[] prevPath2 = { 0, 0 };
	
		for (i = 1; i <= massSize; i++) {
	
			Peak peak = peaks.get(i - 1);
			currentMass = peak.getMass();
			
			currentIntensity = peak.getIntensity();
			currentMassBin = EmitionProb.getBinForMass(currentMass,
					precursorMass); // [self getMassBin:[ [tanMasses
			// objectAtIndex:i-1] floatValue]
			// peptideMass:Pmass];
			currentIntensityBin = EmitionProb.getBinForIntensity(
					currentIntensity, sortedIntensity);
	
			ArrayList<IonMatch> matches = matchCurrentMass(currentMass,
					prevPath, prevPath2, fragsMasses, fragsMatches);
			for (l = 0; l < matches.size(); l++) {
				kmax = 0;
				IonMatch currentMatch = matches.get(l);
				currentIonType = currentMatch.getIonName();
				currentIonIndex = currentMatch.getIonIndex();
				NTerCleaveProb = currentMatch.getNTerCleavage();
				CTerCleaveProb = currentMatch.getCTerCleavage();
	
				maxprod = viterbi[kmax][i - 1] + logA[kmax][currentIonType]; // [viterbi
				// valueAtRow:i-1
				// column:kmax]
				// +[A_ij
				// valueAtRow:kmax
				// column:l
				// ];
				for (k = 1; k < numberOfState; k++) {
					prod = viterbi[k][i - 1] + logA[k][currentIonType];
					if (prod > maxprod) {
						kmax = k;
						maxprod = prod;
					}
				}
				currentMassEmitProb = logEM[currentIonType][currentMassBin]; // [E_M
				// valueAtRow:currentMassBin
				// column:l];
				currentIntensityEmitProb = logEI[currentIonType][currentIntensityBin];
				thisValue = currentMassEmitProb + currentIntensityEmitProb
						+ NTerCleaveProb + CTerCleaveProb + (currentIntensity/totalIntensity);
	
				currentValue = maxprod + thisValue + observedProbArray[currentIonType];
				viterbi[currentIonType][i] = currentValue;
	
				traceBack[i][currentIonType] = kmax;
				currentPath[0] = currentIonType;
				currentPath[1] = currentIonIndex;
			}
			prevPath2[0] = prevPath[0];
			prevPath2[1] = prevPath[1];
			prevPath[0] = currentPath[0];
			prevPath[1] = currentPath[1];
	
		}
	
		kmax = 0;
		maxprod = viterbi[kmax][massSize]; // [ viterbi valueAtRow:massCount
		// column:kmax] ; //v[L][kmax];
	
		for (k = 1; k < numberOfState; k++) {
			prod = viterbi[k][massSize];
			if (prod > maxprod) {
				kmax = k;
				maxprod = prod;
			}
		}
	
		path[massSize - 1] = kmax;
		for (i = massSize - 2; i >= 0; i--) {
			kk = path[i + 1];
			path[i] = traceBack[i + 2][kk];
	
		}
	
		maxprod = (maxprod / 1000);
		maxprod = (Math.exp(maxprod));
		maxprod = maxprod * massSize;
		double weightFactor = calculateWeightFactor(path, massSize, acidSequence);
		maxprod = maxprod * weightFactor;	
		return maxprod;
	}

	/*****************************************************************************************************/
	public static double calculateWeightFactor(int[] path, int size,
			String sequence) {

		int AIonCount = 0, B17IonCount = 0, B18IonCount = 0, Y17IonCount = 0, Y18IonCount = 0, internalIonCount = 0, internal17Count = 0, internal18Count = 0, immIonCount = 0;
		int totalBIon, totalYIon, totalInternalIon;
		double matchPeakFreq = 0.0;
		int BIonCount = 0;
		int YIonCount = 0;

		for (int i = 0; i < size; i++) {
			if (path[i] == Defines.B_ION + 1)
				BIonCount++;
			if (path[i] == Defines.Y_ION + 1)
				YIonCount++;
			if (path[i] == Defines.A_ION + 1)
				AIonCount++;
			if (path[i] == Defines.B17_ION + 1)
				B17IonCount++;
			if (path[i] == Defines.Y17_ION + 1)
				Y17IonCount++;
			if (path[i] == Defines.B18_ION + 1)
				B18IonCount++;
			if (path[i] == Defines.Y18_ION + 1)
				Y18IonCount++;
			if (path[i] == Defines.IMM_ION + 1)
				immIonCount++;
			if (path[i] == Defines.INTERNAL_ION + 1)
				internalIonCount++;
			if (path[i] == Defines.INTERNAL17 + 1)
				internal17Count++;
			if (path[i] == Defines.INTERNAL18 + 1)
				internal18Count++;
		}
		int length = sequence.length();
		totalBIon = length - 1;
		totalYIon = length - 1;
		totalInternalIon = calculateInternalCount(sequence);
		if (totalBIon == 0)
			totalBIon = 1;
		if (totalYIon == 0)
			totalYIon = 1;
		if (totalInternalIon == 0)
			totalInternalIon = 1;
		matchPeakFreq = ((double) BIonCount / totalBIon)
				* observedProbArray[Defines.B_ION]
				+ ((double) YIonCount / totalYIon)
				* observedProbArray[Defines.Y_ION]
				+ ((double) internalIonCount / totalInternalIon)
				* observedProbArray[Defines.INTERNAL_ION]
				+ ((double) B17IonCount / totalBIon)
				* observedProbArray[Defines.B17_ION]
				+ ((double) Y17IonCount / totalYIon)
				* observedProbArray[Defines.Y17_ION]
				+ ((double) internal17Count / totalInternalIon)
				* observedProbArray[Defines.INTERNAL17]
				+ ((double) B18IonCount / totalBIon)
				* observedProbArray[Defines.B18_ION]
				+ ((double) Y18IonCount / totalYIon)
				* observedProbArray[Defines.Y18_ION]
				+ ((double) internal18Count / totalInternalIon)
				* observedProbArray[Defines.INTERNAL18]
				+ ((double) AIonCount / totalBIon)
				* observedProbArray[Defines.A_ION];

		return 5 * matchPeakFreq*length/4;
	}

	/**************************************************************************************/
	public static void compareInternalFrags(Vector ionTypes, Vector NTerProbs,
			Vector CTerProbs, double currentMass, String sequence,
			int[] prevPath, int[] prevPath2, int Type) {
		int seqLen = sequence.length(), i, j;
		char prevAcid, currentAcid, prevCTerAcid, nextNTerAcid;
		int pAcid, cAcid, pCAcid, nNAcid;
		double hydrogenMass = (Type == Defines.MONOISOTOPIC) ? Defines.HYDROGEN_MONO
				: Defines.HYDROGEN_AVE, oxygenMass = (Type == Defines.MONOISOTOPIC) ? Defines.OXYGEN_MONO
				: Defines.OXYGEN_AVE, carbonMass = (Type == Defines.MONOISOTOPIC) ? Defines.CARBON_MONO
				: Defines.CARBON_AVE, nitrogenMass = (Type == Defines.MONOISOTOPIC) ? Defines.NITROGEN_MONO
				: Defines.NITROGEN_AVE;
		double prevAcidMass, currentAcidMass, currentFrag = 0.0;
		double frag17, frag18;
		double NTerCleavage, CTerCleavage;
		int currentIndex = 0;
		for (i = 1; i < seqLen - 2; i++) {
			prevCTerAcid = sequence.charAt(i - 1);
			prevAcid = sequence.charAt(i);
			pCAcid = CleavageProb.getIntValue(prevCTerAcid);
			pAcid = CleavageProb.getIntValue(prevAcid);
			prevAcidMass = HMM_SCORE.AAMasses[prevAcid];
			currentFrag = prevAcidMass + hydrogenMass;
			for (j = i + 1; j < seqLen - 1; j++) {
				currentIndex += i;
				currentAcid = sequence.charAt(j);
				nextNTerAcid = sequence.charAt(j + 1);
				cAcid = CleavageProb.getIntValue(currentAcid);
				nNAcid = CleavageProb.getIntValue(nextNTerAcid);
				currentAcidMass = HMM_SCORE.AAMasses[currentAcid]; // [
				// [residueMasses
				// objectForKey:currentAcid]
				// floatValue];
				currentFrag += currentAcidMass;
				frag17 = currentFrag - nitrogenMass - (3 * hydrogenMass);
				frag18 = currentFrag - oxygenMass - (2 * hydrogenMass);
				if (Math.abs(currentMass - currentFrag) < Properties.fragmentTolerance) {
					if (!(prevPath[0] == Defines.INTERNAL_ION + 1 & prevPath[1] == i)
							& !(prevPath2[0] == Defines.INTERNAL_ION + 1 & prevPath2[1] == i)) {
						NTerCleavage = (logNTerCleavageProb[pAcid] + logNTerCleavageProb[nNAcid]) / 2;
						;
						CTerCleavage = (logCTerCleavageProb[pCAcid] + logCTerCleavageProb[cAcid]) / 2;
						;
						int[] ionName = { Defines.INTERNAL_ION + 1,
								currentIndex };
						ionTypes.add(ionName);
						NTerProbs.add(NTerCleavage);
						CTerProbs.add(CTerCleavage);
					}
				}

				if (Math.abs(currentMass - frag17) < Properties.fragmentTolerance) {
					if (!(prevPath[0] == Defines.INTERNAL17 + 1 & prevPath[1] == i)
							& !(prevPath2[0] == Defines.INTERNAL17 + 1 & prevPath2[1] == i)) {
						NTerCleavage = (logNTerCleavageProb[pAcid] + logNTerCleavageProb[nNAcid]) / 2;
						;
						CTerCleavage = (logCTerCleavageProb[pCAcid] + logCTerCleavageProb[cAcid]) / 2;
						;
						int[] ionName = { Defines.INTERNAL17 + 1, currentIndex };
						ionTypes.add(ionName);
						NTerProbs.add(NTerCleavage);
						CTerProbs.add(CTerCleavage);
					}
				}

				if (Math.abs(currentMass - frag18) < Properties.fragmentTolerance) {
					if (!(prevPath[0] == Defines.INTERNAL18 + 1 & prevPath[1] == i)
							& !(prevPath2[0] == Defines.INTERNAL18 + 1 & prevPath2[1] == i)) {
						NTerCleavage = (logNTerCleavageProb[pAcid] + logNTerCleavageProb[nNAcid]) / 2;
						;
						CTerCleavage = (logCTerCleavageProb[pCAcid] + logCTerCleavageProb[cAcid]) / 2;
						;
						int[] ionName = { Defines.INTERNAL18 + 1, currentIndex };
						ionTypes.add(ionName);
						NTerProbs.add(NTerCleavage);
						CTerProbs.add(CTerCleavage);
					}
				}

			}
		}

	}
	
	
    
    /*************************************************************************************************************/
    public static double getTotalIntensity (ArrayList<Peak> peaks) {
       
	     double totalIntensity=0.0;       
         for (Peak peak : peaks) {
	          totalIntensity += peak.getIntensity();
        }
        return totalIntensity/100;
   }   
	  
	

	/*****************************************************************************************************/
	public static int[] countIonType(String sequence) {

		int totalBIon, totalYIon, totalInternalIon;
		int[] ionCount = new int[numberOfState - 1];

		totalBIon = sequence.length() - 1;
		totalInternalIon = calculateInternalCount(sequence);
		ionCount[Defines.B_ION] = totalBIon;
		ionCount[Defines.Y_ION] = totalBIon;
		ionCount[Defines.B17_ION] = totalBIon;
		ionCount[Defines.B18_ION] = totalBIon;
		ionCount[Defines.Y17_ION] = totalBIon;
		ionCount[Defines.Y18_ION] = totalBIon;
		ionCount[Defines.A_ION] = totalBIon;
		ionCount[Defines.IMM_ION] = totalBIon + 1;
		ionCount[Defines.INTERNAL_ION] = totalInternalIon;
		ionCount[Defines.INTERNAL17] = totalInternalIon;
		ionCount[Defines.INTERNAL18] = totalInternalIon;
		return ionCount;

	}

	/*********************************************************************************************************/
	public static Hashtable getHash(float mass, Hashtable allTypeIons,
			String sequence) {

		int seqLen = sequence.length();
		float peptideMass;
		char Aacid, AnextAcid;
		double acidMass;
		int acid, nextAcid;
		Hashtable matchHash = new Hashtable();
		double bIon = 0.0, prevBIon = 0.0, yIon, aIon, immIon, b17Ion, b18Ion, y17Ion, y18Ion;
		double ionMass;
		Vector ionTypes = new Vector();
		Vector NTerCleaveProbs = new Vector();
		Vector CTerCleaveProbs = new Vector();
		Vector yIons = (Vector) allTypeIons.get("Y_IONS");
		Vector bIons = (Vector) allTypeIons.get("B_IONS");
		Vector aIons = (Vector) allTypeIons.get("A_IONS");
		Vector y17Ions = (Vector) allTypeIons.get("Y17_IONS");
		Vector y18Ions = (Vector) allTypeIons.get("Y18_IONS");
		Vector b17Ions = (Vector) allTypeIons.get("B17_IONS");
		Vector b18Ions = (Vector) allTypeIons.get("B18_IONS");
		Vector internalIons = (Vector) allTypeIons.get("INTERNAL");
		Vector int17Ions = (Vector) allTypeIons.get("INTERNAL17");
		Vector int18Ions = (Vector) allTypeIons.get("INTERNAL18");
		Vector immIons = (Vector) allTypeIons.get("IMM_IONS");
		int i = 0;
		int j;
		boolean matched = false;
		int totalCount = yIons.size();
		while (i < totalCount) {
			Aacid = sequence.charAt(i);
			acid = CleavageProb.getIntValue(Aacid);
			AnextAcid = sequence.charAt(i + 1);
			nextAcid = CleavageProb.getIntValue(AnextAcid);
			j = 0;
			while (j < bIons.size() & !matched) {
				bIon = (Double) bIons.get(j);
				if (Math.abs(mass - bIon) < Properties.fragmentTolerance) {
					ionTypes.add(Defines.B_ION);
					bIons.remove(j);
					NTerCleaveProbs.add(logNTerCleavageProb[nextAcid]);
					CTerCleaveProbs.add(logCTerCleavageProb[acid]);
				}
			}

		}

		if (ionTypes.size() < 1) { // then the current mass did not match to any
			// theoretical ion and put that as noise
			ionTypes.add(Defines.J_ION);
			NTerCleaveProbs.add(Math.log(ZERO_FINAL));
			CTerCleaveProbs.add(Math.log(ZERO_FINAL));
		}
		matchHash.put("IonTypes", ionTypes);
		matchHash.put("NTerCleaveProb", NTerCleaveProbs);
		matchHash.put("CTerCleaveProb", CTerCleaveProbs);
		return matchHash;
	}

	public static int getIntForIonName(String ionName) {
		ionName = ionName.trim();
		int intIon = 1;
		if (ionName.equals("B_ION"))
			intIon = Defines.B_ION;
		else if (ionName.equals("Y_ION"))
			intIon = Defines.Y_ION;
		else if (ionName.equals("IMM_ION"))
			intIon = Defines.IMM_ION;
		else if (ionName.equals("A_ION"))
			intIon = Defines.A_ION;
		else if (ionName.equals("J_ION"))
			intIon = Defines.J_ION;
		else if (ionName.equals("B17_ION"))
			intIon = Defines.B17_ION;
		else if (ionName.equals("Y17_ION"))
			intIon = Defines.Y17_ION;
		else if (ionName.equals("B18_ION"))
			intIon = Defines.B18_ION;
		else if (ionName.equals("Y18_ION"))
			intIon = Defines.Y18_ION;
		else if (ionName.equals("INTERNAL_ION"))
			intIon = Defines.INTERNAL_ION;
		else if (ionName.equals("INTERNAL17"))
			intIon = Defines.INTERNAL17;
		else if (ionName.equals("INTERNAL18"))
			intIon = Defines.INTERNAL18;
		return intIon;
	}

	public static int getIntValue(char acid) {
		if (acid == 'A')
			return Defines.A;
		if (acid == 'L')
			return Defines.L;
		if (acid == 'G')
			return Defines.G;
		if (acid == 'H')
			return Defines.H;
		if (acid == 'M')
			return Defines.M;
		if (acid == 'N')
			return Defines.N;
		if (acid == 'P')
			return Defines.P;
		if (acid == 'D')
			return Defines.D;
		if (acid == 'F')
			return Defines.F;
		if (acid == 'I')
			return Defines.I;
		if (acid == 'K')
			return Defines.K;
		if (acid == 'S')
			return Defines.S;
		if (acid == 'T')
			return Defines.T;
		if (acid == 'W')
			return Defines.W;
		if (acid == 'Y')
			return Defines.Y;
		if (acid == 'V')
			return Defines.V;
		if (acid == 'R')
			return Defines.R;
		if (acid == 'E')
			return Defines.E;
		if (acid == 'Q')
			return Defines.Q;
		if (acid == 'C')
			return Defines.C;
		else
			return Defines.A;
	}

	/*********************************************************************************************************/
	public static Hashtable getIonTypes(double currentMass, int[] prevPath,
			int[] prevPath2, String sequence, int Type) {

		int seqLen = sequence.length();
		float peptideMass;
		char Aacid, AnextAcid;
		double acidMass;
		int acid, nextAcid;
		Hashtable matchHash = new Hashtable();
		double bIon = 0.0, prevBIon = 0.0, yIon, aIon, immIon, b17Ion, b18Ion, y17Ion, y18Ion;
		Vector ionTypes = new Vector();
		Vector ionNumbers = new Vector();
		Vector NTerCleaveProbs = new Vector();
		Vector CTerCleaveProbs = new Vector();
		double hydrogenMass = (Type == Defines.MONOISOTOPIC) ? Defines.HYDROGEN_MONO
				: Defines.HYDROGEN_AVE, oxygenMass = (Type == Defines.MONOISOTOPIC) ? Defines.OXYGEN_MONO
				: Defines.OXYGEN_AVE, carbonMass = (Type == Defines.MONOISOTOPIC) ? Defines.CARBON_MONO
				: Defines.CARBON_AVE, nitrogenMass = (Type == Defines.MONOISOTOPIC) ? Defines.NITROGEN_MONO
				: Defines.NITROGEN_AVE;
		peptideMass = TheoreticalFragment.calculatePeptideMass(sequence);
		for (int i = 0; i < seqLen; i++) {
			Aacid = sequence.charAt(i);
			acid = CleavageProb.getIntValue(Aacid);
			acidMass = HMM_SCORE.AAMasses[Aacid];
			if (i < seqLen - 1) {
				AnextAcid = sequence.charAt(i + 1);
				nextAcid = CleavageProb.getIntValue(AnextAcid);
				if (i == 0)
					bIon = acidMass + hydrogenMass;
				else
					bIon = prevBIon + acidMass;
				if (Math.abs(currentMass - bIon) < Properties.fragmentTolerance) {
					if (!(prevPath[0] == Defines.B_ION + 1 & prevPath[1] == i)
							& !(prevPath2[0] == Defines.B_ION + 1 & prevPath2[1] == i)) {
						int[] ionName = { Defines.B_ION + 1, i };
						ionTypes.add(ionName);
						NTerCleaveProbs.add(logNTerCleavageProb[nextAcid]);
						CTerCleaveProbs.add(logCTerCleavageProb[acid]);
					}
				}
				yIon = peptideMass - bIon + 2 * hydrogenMass;

				if (Math.abs(currentMass - yIon) < Properties.fragmentTolerance) {
					if (!(prevPath[0] == Defines.Y_ION + 1 & prevPath[1] == i)
							& !(prevPath2[0] == Defines.Y_ION + 1 & prevPath2[1] == i)) {
						int[] ionName = { Defines.Y_ION + 1, i };
						ionTypes.add(ionName);
						NTerCleaveProbs.add(logNTerCleavageProb[nextAcid]);
						CTerCleaveProbs.add(logCTerCleavageProb[acid]);
					}
				}
				b17Ion = bIon - nitrogenMass - (3 * hydrogenMass);
				if (Math.abs(currentMass - b17Ion) < Properties.fragmentTolerance) {
					if (!(prevPath[0] == Defines.B17_ION + 1 & prevPath[1] == i)
							& !(prevPath2[0] == Defines.B17_ION + 1 & prevPath2[1] == i)) {
						int[] ionName = { Defines.B17_ION + 1, i };
						ionTypes.add(ionName);
						NTerCleaveProbs.add(logNTerCleavageProb[nextAcid]);
						CTerCleaveProbs.add(logCTerCleavageProb[acid]);
					}
				}
				b18Ion = bIon - oxygenMass - (2 * hydrogenMass);
				if (Math.abs(currentMass - b18Ion) < Properties.fragmentTolerance) {
					if (!(prevPath[0] == Defines.B18_ION + 1 & prevPath[1] == i)
							& !(prevPath2[0] == Defines.B18_ION + 1 & prevPath2[1] == i)) {
						int[] ionName = { Defines.B18_ION + 1, i };
						ionTypes.add(ionName);
						NTerCleaveProbs.add(logNTerCleavageProb[nextAcid]);
						CTerCleaveProbs.add(logCTerCleavageProb[acid]);
					}
				}
				y18Ion = yIon - oxygenMass - (2 * hydrogenMass);
				if (Math.abs(currentMass - y18Ion) < Properties.fragmentTolerance) {
					if (!(prevPath[0] == Defines.Y18_ION + 1 & prevPath[1] == i)
							& !(prevPath2[0] == Defines.Y18_ION + 1 & prevPath2[1] == i)) {
						int[] ionName = { Defines.Y18_ION + 1, i };
						ionTypes.add(ionName);
						NTerCleaveProbs.add(logNTerCleavageProb[nextAcid]);
						CTerCleaveProbs.add(logCTerCleavageProb[acid]);
					}
				}
				y17Ion = yIon - nitrogenMass - (3 * hydrogenMass);
				if (Math.abs(currentMass - y17Ion) < Properties.fragmentTolerance) {
					if (!(prevPath[0] == Defines.Y17_ION + 1 & prevPath[1] == i)
							& !(prevPath2[0] == Defines.Y17_ION + 1 & prevPath2[1] == i)) {
						int[] ionName = { Defines.Y17_ION + 1, i };
						ionTypes.add(ionName);
						NTerCleaveProbs.add(logNTerCleavageProb[nextAcid]);
						CTerCleaveProbs.add(logCTerCleavageProb[acid]);
					}
				}

				aIon = bIon - carbonMass - oxygenMass;
				if (Math.abs(currentMass - aIon) < Properties.fragmentTolerance) {
					if (!(prevPath[0] == Defines.A_ION + 1 & prevPath[1] == i)
							& !(prevPath2[0] == Defines.A_ION + 1 & prevPath2[1] == i)) {
						int[] ionName = { Defines.A_ION + 1, i };
						ionTypes.add(ionName);
						NTerCleaveProbs.add(logNTerCleavageProb[nextAcid]);
						CTerCleaveProbs.add(logCTerCleavageProb[acid]);
					}
				}
			}
			immIon = acidMass - carbonMass - oxygenMass + hydrogenMass;
			if (Math.abs(currentMass - immIon) < Properties.fragmentTolerance) {
				if (!(prevPath[0] == Defines.IMM_ION + 1 & prevPath[1] == i)
						& !(prevPath2[0] == Defines.IMM_ION + 1 & prevPath2[1] == i)) {
					int[] ionName = { Defines.IMM_ION + 1, i };
					ionTypes.add(ionName);
					NTerCleaveProbs.add(logNTerCleavageProb[acid]);
					CTerCleaveProbs.add(logCTerCleavageProb[acid]);
				}
			}

			prevBIon = bIon;

		}
		compareInternalFrags(ionTypes, NTerCleaveProbs, CTerCleaveProbs,
				currentMass, sequence, prevPath, prevPath2, Type);

		if (ionTypes.size() < 1) { // then the current mass did not match to any
			// theoretical ion and put that as noise
			int[] ionName = { Defines.J_ION + 1, 0 };
			ionTypes.add(ionName);
			NTerCleaveProbs.add(Math.log(ZERO_FINAL));
			CTerCleaveProbs.add(Math.log(ZERO_FINAL));
		}
		matchHash.put("IonTypes", ionTypes);
		matchHash.put("NTerCleaveProb", NTerCleaveProbs);
		matchHash.put("CTerCleaveProb", CTerCleaveProbs);

		return matchHash;
	}

	/******************************************************************************************************
	 * Read Probabilities array from file Input: Name of the probability file
	 ********************************************************************************************************/
	private static float[][] getTranProbabilities(String fileName, File dir)
			throws IOException {
		float[][] tranMatrix = new float[Defines.numberOfIons][Defines.numberOfIons];
		int firstIon, secondIon;
		String line = null;
		String ionName1, ionName2 = null;
		float prob;
		String[] lineStr;
		File probFile = new File(dir, fileName);
		if (!probFile.exists()) {
			System.exit(0);
		}
		BufferedReader reader = new BufferedReader(new FileReader(probFile));
		while ((line = reader.readLine()) != null) {
			lineStr = line.split("\t");
			ionName1 = lineStr[0];
			ionName2 = lineStr[1];
			prob = (Float.valueOf(lineStr[2])).floatValue();
			firstIon = getIntForIonName(ionName1);
			secondIon = getIntForIonName(ionName2);
			tranMatrix[firstIon][secondIon] = prob;
		}
		return tranMatrix;

	}

	/******************************************************************************/
	public static void HmmSetUp() {
		File paramDir = Properties.HMMScoreParametersFile;
		float[][] transitionMat = new float[Defines.numberOfIons][Defines.numberOfIons];
		float massEmitMat, intensityEmitMat;
		File[] paramFiles = paramDir.listFiles();
		int fileCount = paramFiles.length;
		try {
			for (int i = 0; i < fileCount; i++) {
				File currentFile = paramFiles[i];
				String currentName = currentFile.getName();
				if (currentName.contentEquals("HMMProbFile")) {
					setProbabilities(currentName, paramDir);
				}
				if (currentName.contentEquals("HMMObsProbFile")) {
					setObsProbabilities(currentName, paramDir);
				}
				if (currentName.contentEquals("HMMTranProbFile")) {
					transitionMat = getTranProbabilities(currentName, paramDir);
				}
				if (currentName.contentEquals("HMMMassEmitProbFile")) {
					setMassEmitProb(currentName, paramDir);
				}
				if (currentName.contentEquals("HMMIntensityEmitProbFile")) {
					setIntensityEmitProb(currentName, paramDir);
				}
				if (currentName.contentEquals("HMMCleavageProbFile")) {
					setCleavageProb(currentName, paramDir);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		setupObservedProb(observedProbArray);
		setupTransitionMatrix(transitionMat, probArray);
		calculateAve();
		Parameter.writeln("HMM Setup done!!!!!!\n");
	}

	/****************************************************************************************/
	public static ArrayList<IonMatch> matchCurrentMass(double mass,
			int[] prevPath, int[] prevPath2, Vector<Double> fragsMasses,
			ArrayList<IonMatch> fragsMatches) {

		IonMatch currentMatch;
		ArrayList<IonMatch> matches = new ArrayList<IonMatch>();
		matches.clear();
		int i, massCount = fragsMasses.size();
		for (i = 0; i < massCount; i++) {
			double currentMass = fragsMasses.get(i);
			if ((mass - currentMass) > -1.0) {
				if (Math.abs(mass - currentMass) < Properties.fragmentTolerance) {
					currentMatch = fragsMatches.get(i);
					int ionName = currentMatch.getIonName();
					int ionIndex = currentMatch.getIonIndex();
					if (!(prevPath[0] == ionName & prevPath[1] == ionIndex)
							& !(prevPath2[0] == ionName & prevPath2[1] == ionIndex)) {
						matches.add(currentMatch);

					}
				}
			}
		}

		if (matches.size() < 1) {
			IonMatch junkIon = new IonMatch(Defines.J_ION + 1, 0, Math
					.log(ZERO_FINAL), Math.log(ZERO_FINAL));
			matches.add(junkIon);
		}
		return matches;
	}

	/*********************************************************************************************************/

	public void generateFragmentsForSequence(String sequence, int Type) {

		double hydrogenMass = (Type == Defines.MONOISOTOPIC) ? Defines.HYDROGEN_MONO
				: Defines.HYDROGEN_AVE, oxygenMass = (Type == Defines.MONOISOTOPIC) ? Defines.OXYGEN_MONO
				: Defines.OXYGEN_AVE, carbonMass = (Type == Defines.MONOISOTOPIC) ? Defines.CARBON_MONO
				: Defines.CARBON_AVE, nitrogenMass = (Type == Defines.MONOISOTOPIC) ? Defines.NITROGEN_MONO
				: Defines.NITROGEN_AVE, waterMass = (Type == Defines.MONOISOTOPIC) ? Defines.WATER_MONO
				: Defines.WATER_AVERAGE, amoniaMass = (Type == Defines.MONOISOTOPIC) ? Defines.AMONIA_MONO
				: Defines.AMONIA_AVERAGE;

		int seqLen = sequence.length();
		int i;
		char Aacid, AnextAcid = 'a';
		double acidMass;
		double bIon = 0.0, prevBIon = 0.0, yIon, aIon, immIon, b17Ion, b18Ion, y17Ion, y18Ion;
		fragsMasses.clear();
		fragsMatches.clear();
		IonMatch currentMatch;
		double peptideMass = calculatePeptideMass(sequence);
		for (i = 0; i < seqLen; i++) {
			Aacid = sequence.charAt(i);
			acidMass = AminoAcids.getWeightMono(Aacid);
			if (i < seqLen - 1) {
				AnextAcid = sequence.charAt(i + 1);
				if (i == 0)
					bIon = acidMass + hydrogenMass;
				else
					bIon = prevBIon + acidMass;
				fragsMasses.add(bIon);
				currentMatch = new IonMatch(Defines.B_ION + 1, i + 1,
						logNTerCleavageProb[AnextAcid],
						logCTerCleavageProb[Aacid]);
				fragsMatches.add(currentMatch);

				int name2 = currentMatch.getIonName();
				int index2 = currentMatch.getIonIndex();
				double cTer2 = currentMatch.getCTerCleavage();
				double nTer2 = currentMatch.getNTerCleavage();

				yIon = peptideMass - bIon + 2 * hydrogenMass;

				fragsMasses.add(yIon);
				currentMatch = new IonMatch(Defines.Y_ION + 1, i + 1,
						logNTerCleavageProb[AnextAcid],
						logCTerCleavageProb[Aacid]);
				fragsMatches.add(currentMatch);

				b17Ion = bIon - amoniaMass;

				fragsMasses.add(b17Ion);
				currentMatch = new IonMatch(Defines.B17_ION + 1, i + 1,
						logNTerCleavageProb[AnextAcid],
						logCTerCleavageProb[Aacid]);
				fragsMatches.add(currentMatch);

				b18Ion = bIon - waterMass;

				fragsMasses.add(b18Ion);
				currentMatch = new IonMatch(Defines.B18_ION + 1, i + 1,
						logNTerCleavageProb[AnextAcid],
						logCTerCleavageProb[Aacid]);
				fragsMatches.add(currentMatch);

				y18Ion = yIon - waterMass;

				fragsMasses.add(y18Ion);
				currentMatch = new IonMatch(Defines.Y18_ION + 1, i + 1,
						logNTerCleavageProb[AnextAcid],
						logCTerCleavageProb[Aacid]);
				fragsMatches.add(currentMatch);

				y17Ion = yIon - amoniaMass;

				fragsMasses.add(y17Ion);
				currentMatch = new IonMatch(Defines.Y17_ION + 1, i + 1,
						logNTerCleavageProb[AnextAcid],
						logCTerCleavageProb[Aacid]);
				fragsMatches.add(currentMatch);

				aIon = bIon - carbonMass - oxygenMass;

				fragsMasses.add(aIon);
				currentMatch = new IonMatch(Defines.A_ION + 1, i + 1,
						logNTerCleavageProb[AnextAcid],
						logCTerCleavageProb[Aacid]);
				fragsMatches.add(currentMatch);

			}
			immIon = acidMass - carbonMass - waterMass;

			fragsMasses.add(immIon);
			currentMatch = new IonMatch(Defines.IMM_ION + 1, i + 1,
					logNTerCleavageProb[AnextAcid], logCTerCleavageProb[Aacid]);
			fragsMatches.add(currentMatch);

			prevBIon = bIon;

		}

		/*
		 * Print the mass and ionTypes for this sequence
		 * 
		 * for (IonMatch match:fragsMatches) { int name1 = match.getIonName();
		 * int index1 = match.getIonIndex(); double cTer1 =
		 * match.getCTerCleavage(); double nTer1 = match.getNTerCleavage(); }
		 */
		calculateInternalFrags(sequence, 0); // sequence:sequence type:Type];

		sortFragments(fragsMasses, fragsMatches);
	}

	/**************************************************************************************/
	private void calculateInternalFrags(String sequence, int Type) {
		int seqLen = sequence.length(), i, j;

		char prevAcid, currentAcid, prevCTerAcid, nextNTerAcid;
		Vector<Double> internalFrags = new Vector<Double>();
		double hydrogenMass = (Type == Defines.MONOISOTOPIC) ? Defines.HYDROGEN_MONO
				: Defines.HYDROGEN_AVE, oxygenMass = (Type == Defines.MONOISOTOPIC) ? Defines.OXYGEN_MONO
				: Defines.OXYGEN_AVE, nitrogenMass = (Type == Defines.MONOISOTOPIC) ? Defines.NITROGEN_MONO
				: Defines.NITROGEN_AVE;

		double prevAcidMass, currentAcidMass, currentFrag = 0.0;
		double frag17, frag18;
		double NTerCleavage, CTerCleavage;
		int currentIndex = 0;
		IonMatch currentMatch;
		for (i = 1; i < seqLen - 2; i++) {
			prevCTerAcid = sequence.charAt(i - 1);
			prevAcid = sequence.charAt(i);
			prevAcidMass = AminoAcids.getWeightMono(prevAcid);
			currentFrag = prevAcidMass + hydrogenMass;
			for (j = i + 1; j < seqLen - 1; j++) {
				currentIndex += i;
				currentAcid = sequence.charAt(j);
				nextNTerAcid = sequence.charAt(j + 1);
				currentAcidMass = AminoAcids.getWeightMono(currentAcid); // [
				// [residueMasses
				// objectForKey:currentAcid]
				// floatValue];
				currentFrag += currentAcidMass;
				fragsMasses.add(currentFrag);
				NTerCleavage = (logNTerCleavageProb[(int) prevAcid] + logNTerCleavageProb[(int) nextNTerAcid]) / 2;
				CTerCleavage = (logCTerCleavageProb[(int) prevCTerAcid] + logCTerCleavageProb[(int) currentAcid]) / 2;
				currentMatch = new IonMatch(Defines.INTERNAL_ION + 1,
						currentIndex, NTerCleavage, CTerCleavage);
				fragsMatches.add(currentMatch);
				frag17 = currentFrag - nitrogenMass - (3 * hydrogenMass);
				fragsMasses.add(frag17);
				NTerCleavage = (logNTerCleavageProb[(int) prevAcid] + logNTerCleavageProb[(int) nextNTerAcid]) / 2;
				CTerCleavage = (logCTerCleavageProb[(int) prevCTerAcid] + logCTerCleavageProb[(int) currentAcid]) / 2;
				currentMatch = new IonMatch(Defines.INTERNAL17 + 1,
						currentIndex, NTerCleavage, CTerCleavage);
				fragsMatches.add(currentMatch);

				frag18 = currentFrag - oxygenMass - (2 * hydrogenMass);
				fragsMasses.add(frag18);
				NTerCleavage = (logNTerCleavageProb[(int) prevAcid] + logNTerCleavageProb[(int) nextNTerAcid]) / 2;
				CTerCleavage = (logCTerCleavageProb[(int) prevCTerAcid] + logCTerCleavageProb[(int) currentAcid]) / 2;
				currentMatch = new IonMatch(Defines.INTERNAL18 + 1,
						currentIndex, NTerCleavage, CTerCleavage);
				fragsMatches.add(currentMatch);

			}
		}

	}

	/*******************************************************************************/
	public void sortFragments(Vector<Double> masses, ArrayList<IonMatch> matches) {

		int i, massSize = masses.size();
		double tempMass, currentMass;
		int position;
		IonMatch currentMatch, tempMatch;
		for (i = 1; i < massSize; i++) {
			currentMass = masses.get(i);
			currentMatch = matches.get(i);
			position = i;
			while (position > 0 && (masses.get(position - 1) > currentMass)) {
				tempMass = masses.get(position - 1);
				masses.set(position, tempMass);
				tempMatch = matches.get(position - 1);
				matches.set(position, tempMatch);
				position--;
			}
			masses.set(position, currentMass);
			matches.set(position, currentMatch);
		}
	}

	/******************************************************************************************************
	 * Read Probabilities array from file Input: Name of the probability file
	 ********************************************************************************************************/
	private static void setCleavageProb(String fileName, File dir)
			throws IOException {
		float NTerCleavage, CTerCleavage;
		String line = null;
		String acidName = null;
		float prob;
		String[] lineStr;
		File probFile = new File(dir, fileName);
		if (!probFile.exists()) {
			System.exit(0);
		}
		BufferedReader reader = new BufferedReader(new FileReader(probFile));
		while ((line = reader.readLine()) != null) {
			lineStr = line.split("\t");
			acidName = lineStr[0];
			char acidChar = acidName.charAt(0);
			NTerCleavage = (Float.valueOf(lineStr[1])).floatValue();
			CTerCleavage = (Float.valueOf(lineStr[2])).floatValue();
			if (NTerCleavage < ZERO_FINAL)
				NTerCleavage = ZERO_FINAL;
			if (CTerCleavage < ZERO_FINAL)
				CTerCleavage = ZERO_FINAL;
			logNTerCleavageProb[acidChar] = Math.log(NTerCleavage);
			logCTerCleavageProb[acidChar] = Math.log(CTerCleavage);

		}

	}

	/******************************************************************************************************
	 * Read Probabilities array from file Input: Name of the probability file
	 ********************************************************************************************************/
	private static void setIntensityEmitProb(String fileName, File dir)
			throws IOException {

		int ionType;
		String line = null;
		String ionName = null;
		float prob;
		String[] lineStr;
		File probFile = new File(dir, fileName);
		if (!probFile.exists()) {
			System.exit(0);
		}
		BufferedReader reader = new BufferedReader(new FileReader(probFile));
		while ((line = reader.readLine()) != null) {
			lineStr = line.split("\t");
			ionName = lineStr[0];
			int intIon = getIntForIonName(ionName);
			for (int j = 0; j < Defines.bin_count; j++) {
				logEM[0][j] = Double.NEGATIVE_INFINITY;
				double currentProb = (Float.valueOf(lineStr[j + 1]))
						.floatValue();
				currentProb = Math.log(currentProb);
				logEI[intIon + 1][j] = currentProb;
			}

		}

	}

	/******************************************************************************************************
	 * Read Probabilities array from file Input: Name of the probability file
	 ********************************************************************************************************/
	private static void setMassEmitProb(String fileName, File dir)
			throws IOException {
		int ionType;
		String line = null;
		String ionName = null;
		float prob;
		String[] lineStr;
		File probFile = new File(dir, fileName);
		if (!probFile.exists()) {
			System.exit(0);
		}
		BufferedReader reader = new BufferedReader(new FileReader(probFile));
		while ((line = reader.readLine()) != null) {
			lineStr = line.split("\t");
			ionName = lineStr[0];
			int intIon = getIntForIonName(ionName);
			for (int j = 0; j < Defines.bin_count; j++) {
				logEM[0][j] = Double.NEGATIVE_INFINITY;
				double currentProb = (Float.valueOf(lineStr[j + 1]))
						.floatValue();
				currentProb = Math.log(currentProb);
				logEM[intIon + 1][j] = currentProb;
			}

		}

	}

	/******************************************************************************************************
	 * Read Probabilities array from file Input: Name of the probability file
	 ********************************************************************************************************/
	private static void setObsProbabilities(String fileName, File dir)
			throws IOException {
		String line = null;
		String ionName = null;
		float prob;
		String[] lineStr;
		File probFile = new File(dir, fileName);
		if (!probFile.exists()) {
			System.exit(0);
		}
		BufferedReader reader = new BufferedReader(new FileReader(probFile));
		while ((line = reader.readLine()) != null) {
			lineStr = line.split(" ");
			ionName = lineStr[0];
			prob = (Float.valueOf(lineStr[1])).floatValue();
			int intIon = getIntForIonName(ionName);
			observedProbArray[intIon] = prob;

		}

	}

	/******************************************************************************************************
	 * Read Probabilities array from file Input: Name of the probability file
	 ********************************************************************************************************/
	private static void setProbabilities(String fileName, File dir)
			throws IOException {
		String line = null;
		String ionName = null;
		float prob;
		String[] lineStr;
		File probFile = new File(dir, fileName);
		if (!probFile.exists()) {
			System.exit(0);
		}
		BufferedReader reader = new BufferedReader(new FileReader(probFile));
		while ((line = reader.readLine()) != null) {
			lineStr = line.split(" ");
			ionName = lineStr[0];
			prob = (Float.valueOf(lineStr[1])).floatValue();
			int intIon = getIntForIonName(ionName);
			probArray[intIon] = prob;

		}

	}

	/***********************************************************************************************/
	public static void setupObservedProb(float[] observedProb) {

		int i;
		logObservedProb[0] = Double.NEGATIVE_INFINITY;
		for (i = 1; i < numberOfState; i++) {
			if (observedProb[i - 1] < ZERO_FINAL)
				observedProb[i - 1] = ZERO_FINAL;
			logObservedProb[i] = Math.log(observedProb[i - 1]);
		}
	}

	/**********************************************************************************************/
	public static void setupTransitionMatrix(float[][] tranMatrix, float[] prob) {
		int i, j;
		logA[0][0] = Double.NEGATIVE_INFINITY; // log (0) is -infinity
		for (j = 1; j < numberOfState; j++) {
			if (prob[j - 1] < ZERO_FINAL)
				prob[j - 1] = ZERO_FINAL;
			logA[0][j] = Math.log(prob[j - 1]);
		}
		for (i = 1; i < numberOfState; i++) {
			logA[i][0] = Double.NEGATIVE_INFINITY; // log (0) is -infinity
			for (j = 1; j < numberOfState; j++) {
				if (tranMatrix[i - 1][j - 1] < ZERO_FINAL)
					tranMatrix[i - 1][j - 1] = ZERO_FINAL;
				logA[i][j] = Math.log(tranMatrix[i - 1][j - 1]);

			}
		}
	}

	public static double calculatePeptideMass(String sequence) {
		double acidMass;
		float peptideMass = 0;
		int i, seqLen = sequence.length();
		char acid;
		for (i = 0; i < seqLen; i++) {
			acid = sequence.charAt(i); // [sequence objectAtIndex:i];
			acidMass = AminoAcids.getWeightMono(acid);
			; // [residueMasses objectForKey:acid] ;
			peptideMass += acidMass;
		}
		peptideMass += (Defines.OXYGEN_MONO) + 2 * (Defines.HYDROGEN_MONO);
		return peptideMass;
	}

}
