package HMMScore;
import java.util.Hashtable;
import java.util.Vector;

import Peppy.Properties;

public class CleavageProb {

	private static float[] NTerCleaveProb = new float[Defines.numberOfAminoAcid];
	private static float[] CTerCleaveProb = new float[Defines.numberOfAminoAcid];

	/***********************************************************************************************************************************/
	public static void calculateCleavageProbabilityExcludeTerminal(
			Vector tanMassesArray, Vector acidSequencesArray, int Type) {

		double hydrogenMass = (Type == Defines.MONOISOTOPIC) ? Defines.HYDROGEN_MONO
				: Defines.HYDROGEN_AVE, oxygenMass = (Type == Defines.MONOISOTOPIC) ? Defines.OXYGEN_MONO
				: Defines.OXYGEN_AVE, carbonMass = (Type == Defines.MONOISOTOPIC) ? Defines.CARBON_MONO
				: Defines.CARBON_AVE, nitrogenMass = (Type == Defines.MONOISOTOPIC) ? Defines.NITROGEN_MONO
				: Defines.NITROGEN_AVE;
		int seqLen;
		int i, k, j;
		TheoreticalFragment fragment = new TheoreticalFragment(); // [
																	// [TheoreticalFragment
																	// alloc]
																	// init];
		int massIndex = 0, currentMassCount = 0, tanMassesCount = tanMassesArray
				.size();
		double acidMass, prevBIon, peptideMass;
		char acid, nAcid, firstAcid;
		boolean existsIon = false;
		int[] NTerCleaveCount = new int[Defines.numberOfAminoAcid];
		int[] CTerCleaveCount = new int[Defines.numberOfAminoAcid];
		int[] totalAcidCount = new int[Defines.numberOfAminoAcid];
		for (k = 0; k < Defines.numberOfAminoAcid; k++) {
			totalAcidCount[k] = 0;
			NTerCleaveCount[k] = 0;
			CTerCleaveCount[k] = 0;
		}
		double bIon, yIon, aIon, immIon, b17Ion, b18Ion, y17Ion, y18Ion, a17Ion, a18Ion, cIon, zIon, xIon;
		for (i = 0; i < tanMassesCount; i++) {
			prevBIon = 0.0;
			Hashtable tandemHash = (Hashtable) tanMassesArray.get(i);
			Vector masses = (Vector) tandemHash.get("Masses");
			String acidSequence = (String) acidSequencesArray.get(i);
			seqLen = acidSequence.length();
			peptideMass = fragment.calculatePeptideMass(acidSequence);
			firstAcid = acidSequence.charAt(0);
			for (j = 1; j < seqLen - 2; j++) {
				acid = acidSequence.charAt(j);
				nAcid = acidSequence.charAt(j + 1);
				acidMass = HMM_SCORE.AAMasses[acid]; // [[residueMasses
														// objectForKey:acid]
														// floatValue];
				int currentAcid = getIntValue(acid);
				int nextAcid = getIntValue(nAcid);
				totalAcidCount[currentAcid]++;
				if (j == seqLen - 3)
					totalAcidCount[nextAcid]++;
				if (j == 1) {
					bIon = acidMass + HMM_SCORE.AAMasses[firstAcid]
							+ hydrogenMass;
				} else {
					bIon = prevBIon + acidMass;
				}
				prevBIon = bIon;
				existsIon = checkIon(bIon, masses);
				if (existsIon) {
					CTerCleaveCount[currentAcid]++;
					NTerCleaveCount[nextAcid]++;
				}

				yIon = peptideMass - bIon + 2 * hydrogenMass;
				if (!existsIon) {
					existsIon = checkIon(yIon, masses);
					if (existsIon) {
						CTerCleaveCount[currentAcid]++;
						NTerCleaveCount[nextAcid]++;
					}
				}

				if (!existsIon) {
					b17Ion = bIon - nitrogenMass - (3 * hydrogenMass);
					existsIon = checkIon(b17Ion, masses);
					if (existsIon) {
						CTerCleaveCount[currentAcid]++;
						NTerCleaveCount[nextAcid]++;
					}
				}

				if (!existsIon) {
					b18Ion = bIon - oxygenMass - (2 * hydrogenMass);
					existsIon = checkIon(b18Ion, masses);
					if (existsIon) {
						CTerCleaveCount[currentAcid]++;
						NTerCleaveCount[nextAcid]++;
					}
				}

				if (!existsIon) {
					aIon = bIon - carbonMass - oxygenMass;
					existsIon = checkIon(aIon, masses);
					if (existsIon) {
						CTerCleaveCount[currentAcid]++;
						NTerCleaveCount[nextAcid]++;
					}
				}

				if (!existsIon) {
					a17Ion = bIon - carbonMass - oxygenMass - nitrogenMass
							- (3 * hydrogenMass);
					existsIon = checkIon(a17Ion, masses);
					if (existsIon) {
						CTerCleaveCount[currentAcid]++;
						NTerCleaveCount[nextAcid]++;
					}
				}

				if (!existsIon) {
					a18Ion = bIon - carbonMass - oxygenMass - oxygenMass
							- (2 * hydrogenMass);
					existsIon = checkIon(a18Ion, masses);
					if (existsIon) {
						CTerCleaveCount[currentAcid]++;
						NTerCleaveCount[nextAcid]++;
					}
				}
				if (!existsIon) {
					cIon = bIon + nitrogenMass + 3 * hydrogenMass;
					existsIon = checkIon(cIon, masses);
					if (existsIon) {
						CTerCleaveCount[currentAcid]++;
						NTerCleaveCount[nextAcid]++;
					}
				}

				if (!existsIon) {
					y17Ion = yIon - nitrogenMass - (3 * hydrogenMass);
					existsIon = checkIon(y17Ion, masses);
					if (existsIon) {
						CTerCleaveCount[currentAcid]++;
						NTerCleaveCount[nextAcid]++;
					}
				}

				if (!existsIon) {
					y18Ion = yIon - oxygenMass - (2 * hydrogenMass);
					existsIon = checkIon(y18Ion, masses);
					if (existsIon) {
						CTerCleaveCount[currentAcid]++;
						NTerCleaveCount[nextAcid]++;
					}
				}

				if (!existsIon) {
					zIon = peptideMass - bIon - nitrogenMass + hydrogenMass;
					existsIon = checkIon(zIon, masses);
					if (existsIon) {
						CTerCleaveCount[currentAcid]++;
						NTerCleaveCount[nextAcid]++;
					}
				}
				if (!existsIon) {
					xIon = peptideMass - bIon + carbonMass + oxygenMass;
					existsIon = checkIon(xIon, masses);
					if (existsIon) {
						CTerCleaveCount[currentAcid]++;
						NTerCleaveCount[nextAcid]++;
					}
				}
			}
		}

		for (k = 0; k < Defines.numberOfAminoAcid; k++) {
			NTerCleaveProb[k] = (float) NTerCleaveCount[k] / totalAcidCount[k];
			CTerCleaveProb[k] = (float) CTerCleaveCount[k] / totalAcidCount[k];
			char c = getCharForInt(k);
		}

	}

	public static boolean checkIon(double bIon, Vector masses) {
		int i;
		boolean matched = false;
		int massCount = masses.size();
		for (i = 0; i < massCount; i++) {
			float massValue = (Float) masses.get(i);
			if (Math.abs(massValue - bIon) < Properties.peakDifferenceThreshold) {
				matched = true;
				return matched;
			}
		}
		return matched;
	}

	/***************************************************************************************************************/
	public static float[] CTerCleaveProb() {
		return CTerCleaveProb;
	}

	public static char getCharForInt(int i) {
		if (i == 0)
			return 'A';
		if (i == 1)
			return 'R';
		if (i == 2)
			return 'N';
		if (i == 3)
			return 'D';
		if (i == 4)
			return 'C';
		if (i == 5)
			return 'E';
		if (i == 6)
			return 'Q';
		if (i == 7)
			return 'G';
		if (i == 8)
			return 'H';
		if (i == 9)
			return 'I';
		if (i == 10)
			return 'L';
		if (i == 11)
			return 'K';
		if (i == 12)
			return 'M';
		if (i == 13)
			return 'F';
		if (i == 14)
			return 'P';
		if (i == 15)
			return 'S';
		if (i == 16)
			return 'T';
		if (i == 17)
			return 'W';
		if (i == 18)
			return 'Y';
		if (i == 19)
			return 'V';
		else
			return 'A';
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

	/***************************************************************************************************************/
	public static float[] NTerCleaveProb() {
		return NTerCleaveProb;
	}

	public static void write(String s) {
	}

}
