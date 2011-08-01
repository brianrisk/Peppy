package HMMScore;
import java.util.Hashtable;
import java.util.Vector;

import Peppy.Properties;

public class TheoreticalFragment {

	static Hashtable allTypeIons;

	/***************************************************************************************/
	public static Vector calculateAmoniaLossInternal(Vector internalFrags,
			int Type) {
		int fragCount = internalFrags.size();
		int i;
		Vector retArray = new Vector(); // [NSMutableArray array];
		double hydrogenMass = (Type == Defines.MONOISOTOPIC) ? Defines.HYDROGEN_MONO
				: Defines.HYDROGEN_AVE, oxygenMass = (Type == Defines.MONOISOTOPIC) ? Defines.OXYGEN_MONO
				: Defines.OXYGEN_AVE, carbonMass = (Type == Defines.MONOISOTOPIC) ? Defines.CARBON_MONO
				: Defines.CARBON_AVE, nitrogenMass = (Type == Defines.MONOISOTOPIC) ? Defines.NITROGEN_MONO
				: Defines.NITROGEN_AVE;
		double mass;
		for (i = 0; i < fragCount; i++) {
			mass = (Double) internalFrags.get(i);
			mass = mass - nitrogenMass - (3 * hydrogenMass);
			retArray.add(mass);
		}
		return retArray;
	}

	/**************************************************************************************/
	public static Vector calculateInternalFrags(String sequence, Vector bIons,
			int Type) {
		int seqLen = sequence.length(), i, j;
		char prevAcid, currentAcid;
		boolean isFrag = false;
		double hydrogenMass = (Type == Defines.MONOISOTOPIC) ? Defines.HYDROGEN_MONO
				: Defines.HYDROGEN_AVE;
		double prevAcidMass, currentAcidMass, currentFrag = 0.0;
		Vector internalFrags = new Vector();
		for (i = 1; i < seqLen - 2; i++) {
			prevAcid = sequence.charAt(i);
			prevAcidMass = HMM_SCORE.AAMasses[prevAcid];
			currentFrag = prevAcidMass + hydrogenMass;
			for (j = i + 1; j < seqLen - 1; j++) {
				currentAcid = sequence.charAt(j);
				currentAcidMass = HMM_SCORE.AAMasses[currentAcid]; // [
																	// [residueMasses
																	// objectForKey:currentAcid]
																	// floatValue];
				currentFrag += currentAcidMass;
				isFrag = checkFrag(internalFrags, bIons, currentFrag);
				if (!isFrag)
					internalFrags.add(currentFrag);
			}
		}

		return internalFrags;

	}

	public static float calculatePeptideMass(String sequence) {
		double acidMass;
		float peptideMass = 0;
		int i, seqLen = sequence.length();
		char acid;
		for (i = 0; i < seqLen; i++) {
			acid = sequence.charAt(i); // [sequence objectAtIndex:i];
			acidMass = HMM_SCORE.AAMasses[acid]; // [residueMasses
													// objectForKey:acid] ;
			peptideMass += acidMass;
		}
		peptideMass += (Defines.OXYGEN_MONO) + 2 * (Defines.HYDROGEN_MONO);
		return peptideMass;
	}

	/***************************************************************************************/
	public static Vector calculateWaterLossInternal(Vector internalFrags,
			int Type) {
		int fragCount = internalFrags.size();
		int i;
		Vector retArray = new Vector(); // [NSMutableArray array];
		double hydrogenMass = (Type == Defines.MONOISOTOPIC) ? Defines.HYDROGEN_MONO
				: Defines.HYDROGEN_AVE, oxygenMass = (Type == Defines.MONOISOTOPIC) ? Defines.OXYGEN_MONO
				: Defines.OXYGEN_AVE, carbonMass = (Type == Defines.MONOISOTOPIC) ? Defines.CARBON_MONO
				: Defines.CARBON_AVE, nitrogenMass = (Type == Defines.MONOISOTOPIC) ? Defines.NITROGEN_MONO
				: Defines.NITROGEN_AVE;
		double mass;
		for (i = 0; i < fragCount; i++) {
			mass = (Double) internalFrags.get(i);
			mass = mass - oxygenMass - (2 * hydrogenMass);
			retArray.add(mass);
		}
		return retArray;
	}

	/***************************************************************************************************/
	public static boolean checkFrag(Vector internalFrags, Vector bIons,
			double currentFrag) {
		int i, j;
		boolean isFrag = false;
		for (i = 0; i < internalFrags.size(); i++) {
			double currentMass = (Double) internalFrags.get(i);
			if (Math.abs(currentMass - currentFrag) < Properties.fragmentTolerance) {
				isFrag = true;
				return isFrag;
			}
		}

		for (i = 0; i < bIons.size(); i++) {
			double currentMass = (Double) bIons.get(i);
			if (Math.abs(currentMass - currentFrag) < Properties.fragmentTolerance) {
				isFrag = true;
				return isFrag;
			}
		}

		return isFrag;
	}

	public static Hashtable generateFragmentOfIonType(String sequence, int Type) {

		int seqLen = sequence.length();
		char acid;
		double acidMass;
		double bIon = 0.0, yIon, aIon, immIon, b17Ion, b18Ion, y17Ion, y18Ion, a17Ion, a18Ion, cIon, zIon, xIon;
		double peptide17, peptide18, prevBIon = 0.0, peptideMass;

		Vector bIons = new Vector(), aIons = new Vector(), cIons = new Vector(), zIons = new Vector();
		Vector yIons = new Vector(), b17Ions = new Vector(), b18Ions = new Vector(), y17Ions = new Vector();
		Vector y18Ions = new Vector(), a17Ions = new Vector(), a18Ions = new Vector(), xIons = new Vector();
		Vector immIons = new Vector(), internalIons = new Vector(), internal17Ions = new Vector(), internal18Ions = new Vector();
		Vector unASssined = new Vector();
		double hydrogenMass = (Type == Defines.MONOISOTOPIC) ? Defines.HYDROGEN_MONO
				: Defines.HYDROGEN_AVE, oxygenMass = (Type == Defines.MONOISOTOPIC) ? Defines.OXYGEN_MONO
				: Defines.OXYGEN_AVE, carbonMass = (Type == Defines.MONOISOTOPIC) ? Defines.CARBON_MONO
				: Defines.CARBON_AVE, nitrogenMass = (Type == Defines.MONOISOTOPIC) ? Defines.NITROGEN_MONO
				: Defines.NITROGEN_AVE;
		peptideMass = calculatePeptideMass(sequence);
		for (int i = 0; i < seqLen; i++) {
			acid = sequence.charAt(i);
			acidMass = HMM_SCORE.AAMasses[acid];
			if (i < seqLen - 1) {
				if (i == 0)
					bIon = acidMass + hydrogenMass;
				else
					bIon = prevBIon + acidMass;
				bIons.add(bIon);
				b17Ion = bIon - nitrogenMass - (3 * hydrogenMass);
				b17Ions.add(b17Ion);
				b18Ion = bIon - oxygenMass - (2 * hydrogenMass);
				b18Ions.add(b18Ion);
				yIon = peptideMass - bIon + 2 * hydrogenMass;
				yIons.add(yIon);
				y17Ion = yIon - nitrogenMass - (3 * hydrogenMass);
				y17Ions.add(y17Ion);
				y18Ion = yIon - oxygenMass - (2 * hydrogenMass);
				y18Ions.add(y18Ion);

				aIon = bIon - carbonMass - oxygenMass;
				aIons.add(aIon);

			}
		}
		return allTypeIons;
	}

	public TheoreticalFragment() {
		allTypeIons = new Hashtable();
	}

}
