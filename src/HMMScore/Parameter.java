package HMMScore;
import java.util.Vector;

import Peppy.Properties;

public class Parameter {

	static float[] probArray; // = new float[Defines.numberOfIons];
	static float[] observedProbArray;
	static float[] missProbArray;// = new float[Defines.numberOfIons];
	static float[][] massEmitMat;// = new
									// float[Defines.bin_count][Defines.numberOfIons];
	static float[][] intensityEmitMat;// = new
										// float[Defines.bin_count][Defines.numberOfIons];
	static float[][] transitionMatrix;// = new
										// float[Defines.numberOfIons][Defines.numberOfIons];
	static float[] NTerCleaveProb;// = new float[Defines.numberOfAminoAcid];
	static float[] CTerCleaveProb;// = new float[Defines.numberOfAminoAcid];

	public static float[] getCTerCleavage() {
		return CTerCleaveProb;
	}

	public static float[][] getIntensityEmitMat() {
		return intensityEmitMat;
	}

	public static int getIonTypeCount(Vector bTypeIons, Vector masses) {
		int i, j, bIonCount = 0;
		double typeValue = 0.0;
		for (i = 0; i < masses.size(); i++) {
			float massValue = (Float) masses.get(i);
			for (j = 0; j < bTypeIons.size(); j++) {
				typeValue = (Double) bTypeIons.get(j); // .valueOf(typeStr)).floatValue();
				if (Math.abs(massValue - typeValue) < Properties.fragmentTolerance) {
					bIonCount++;
				}
			}
		}
		return bIonCount;
	}

	public static String getIonTypeForInt(int i) {
		if (i == 0)
			return "B_ION";
		if (i == 1)
			return "Y_ION";
		if (i == 2)
			return "A_ION";
		if (i == 3)
			return "J_ION";
		if (i == 4)
			return "B17_ION";
		if (i == 5)
			return "Y17_ION";
		if (i == 6)
			return "B18_ION";
		if (i == 7)
			return "Y18_ION";
		if (i == 8)
			return "IMM-ION";
		if (i == 9)
			return "INTERNAL_ION";
		if (i == 10)
			return "INTERNAL17";
		if (i == 11)
			return "INTERNAL18";
		if (i == 12)
			return "MISSING";
		if (i == 13)
			return "C_ION";
		if (i == 14)
			return "Z_ION";
		if (i == 15)
			return "X_ION";
		if (i == 16)
			return "A18_ION";
		if (i == 17)
			return "A17_ION";
		else
			return "Y_ION";
	}

	public static float[][] getMassEmitMat() {
		return massEmitMat;
	}

	public static float[] getMissProbabilities() {
		return missProbArray;
	}

	public static float[] getNTerCleavage() {
		return NTerCleaveProb;
	}

	public static float[] getObservedProbabilities() {
		return observedProbArray;
	}

	public static float[] getProbabilities() {
		return probArray;
	}

	public static float[][] getTransitionMat() {
		return transitionMatrix;
	}

	public static void write(String s) {
	}

	public static void writeCleavageProb(float[] prob) {
		int i;
		write("Writing Cleavage Probabilities of each Amino Acid\n");
		for (i = 0; i < Defines.numberOfAminoAcid; i++) {
			char aminoAcid = CleavageProb.getCharForInt(i);
			writeln("Here Amino Acid  and probabilities are " + aminoAcid
					+ " \t" + prob[i]);
		}
	}

	public static void writeln(String s) {
	}

	public static void writeMassEmitMatrix(float[][] massEmitMat) {
		int i, k;
		for (i = Defines.numberOfIons - 1; i >= 0; i--) {
			write("\t" + i);
			for (k = 0; k < Defines.bin_count; k++) {
				write("\t" + massEmitMat[i][k]); // [massEmitMat valueAtRow:i
													// column:k]);
			}
			write("\n");

		}
	}

	public static void writeProbabilities(float[] prob) {
		int i;
		write("Writing Probabilities of each Ion Type\n");
		for (i = 0; i < Defines.numberOfIons; i++) {
			String ionType = getIonTypeForInt(i);
			writeln("Here ion type and probabilities are " + ionType + " \t"
					+ prob[i]);
		}
	}

	public Parameter() {
		observedProbArray = new float[Defines.numberOfIons];
		probArray = new float[Defines.numberOfIons];
		missProbArray = new float[Defines.numberOfIons];
		massEmitMat = new float[Defines.bin_count][Defines.numberOfIons];
		intensityEmitMat = new float[Defines.bin_count][Defines.numberOfIons];
		transitionMatrix = new float[Defines.numberOfIons][Defines.numberOfIons];
		NTerCleaveProb = new float[Defines.numberOfAminoAcid];
		CTerCleaveProb = new float[Defines.numberOfAminoAcid];

	}

}
