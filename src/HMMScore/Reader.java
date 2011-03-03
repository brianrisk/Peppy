package HMMScore;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Hashtable;
import java.util.Vector;

public class Reader {

	public static double calculateMassForCharge(int type, int charge,
			float mOverZ) {
		double mass = 0;
		double protonMass;
		switch (type) {
		case Defines.MONOISOTOPIC:
			protonMass = Defines.HYDROGEN_MONO;
			break;
		case Defines.AVERAGE:
			protonMass = Defines.HYDROGEN_AVE;
			break;
		default:
			protonMass = Defines.HYDROGEN_MONO;
		}

		if (charge > 0)
			mass = charge * (mOverZ - protonMass);
		else if (charge == 0)
			mass = mOverZ;
		return mass;
	}

	/***********************************************************************************************/
	public static Hashtable cleanPeakForIntensity(Hashtable tandemData)
			throws java.io.IOException {
		Hashtable tandemDict = new Hashtable();
		Vector retMasses = new Vector();
		Vector retIntensities = new Vector();
		int maxIntensityPos, massCount;
		int numberToDelete = 0;
		float currentMass, currentIntensity, prevIntensity, maxIntensity, tempIntensity, tempMass;
		Vector masses = (Vector) tandemData.get("Masses");
		massCount = masses.size();
		if (massCount > Defines.TOTAL_PEAK_COUNT) {
			float precursorMass = (Float) tandemData.get("PrecursorMass");
			float precursorIntensity = (Float) tandemData
					.get("PrecursorIntensity");
			int precursorCharge = (Integer) tandemData.get("PrecursorCharge");
			tandemDict.put("PrecursorMass", precursorMass);
			tandemDict.put("PrecursorCharge", precursorCharge);
			tandemDict.put("PrecursorIntensity", precursorIntensity);
			Vector intensities = (Vector) tandemData.get("Intensities");
			for (int i = 0; i < massCount; i++) {
				prevIntensity = (Float) intensities.get(i);
				maxIntensityPos = i;
				for (int j = i + 1; j < massCount; j++) {
					maxIntensity = (Float) intensities.get(maxIntensityPos);
					currentIntensity = (Float) intensities.get(j);
					if (currentIntensity > maxIntensity) {
						maxIntensityPos = j;

					}
				}
				tempIntensity = (Float) intensities.get(maxIntensityPos);
				intensities.setElementAt(intensities.get(i), maxIntensityPos);
				intensities.setElementAt(tempIntensity, i);
				tempMass = (Float) masses.get(maxIntensityPos);
				masses.setElementAt(masses.get(i), maxIntensityPos);
				masses.setElementAt(tempMass, i);

			}
			for (int i = 0; i < Defines.TOTAL_PEAK_COUNT; i++) {
				retMasses.add(masses.get(i));
				retIntensities.add(intensities.get(i));
			}
			tandemDict.put("Masses", retMasses);
			tandemDict.put("Intensities", retIntensities);
			return tandemDict;
		} else
			return tandemData;
	}

	/**************************************************************************************************/

	public static void printTandemData(Hashtable tandemData) {
		Vector masses = (Vector) tandemData.get("Masses");
		Vector intensities = (Vector) tandemData.get("Intensities");
		int i, massCount = masses.size(); // NSNumber *mass;
		float mass, intensity;
		for (i = 0; i < massCount; i++) {
			mass = (Float) masses.get(i);
			intensity = (Float) intensities.get(i);
			writeln("HERE MASS and INTENSITY" + mass + "    " + intensity);
		}
	}

	/**********************************************************************************/
	/*
	 * Read residue's sequence from a file and return as NSArray of NSString
	 * *****
	 * *********************************************************************
	 * ********
	 */

	public static String readSequence(File filename) throws java.io.IOException {
		StringBuffer fileData = new StringBuffer(1000);
		BufferedReader reader = new BufferedReader(new FileReader(filename));
		char[] buf = new char[1024];
		int numRead = 0;
		while ((numRead = reader.read(buf)) != -1) {
			String readData = String.valueOf(buf, 0, numRead);
			fileData.append(readData);
			buf = new char[1024];
		}
		reader.close();
		return fileData.toString();
	}

	/****************************************************************************************************
	 * /* Read residue's masses and name from the file and return as Dictionary
	 * with mass and intensity
	 *****************************************************************************************************/
	public static Hashtable readTandemData(File filename)
			throws java.io.IOException {
		String line = null;
		String massString, intensityStr;
		BufferedReader reader = new BufferedReader(new FileReader(filename));
		String[] lineStr;
		float precursorMZ = 0, precursorIntensity = 0;
		int precursorCharge = 0;
		float precursorMass = 0;
		Hashtable tandemDict = new Hashtable();
		Vector masses = new Vector();
		Vector intensities = new Vector();
		float mOverZ, intensity;

		line = reader.readLine();
		lineStr = line.split(" ");
		if (lineStr.length > 2) {
			precursorMZ = (Float.valueOf(lineStr[0])).floatValue();
			precursorIntensity = (Float.valueOf(lineStr[1])).floatValue();
			precursorCharge = (Integer.valueOf(lineStr[2])).intValue();
			precursorMass = (float) calculateMassForCharge(0, precursorCharge,
					precursorMZ);
		} else if (lineStr.length > 1) {
			precursorMZ = (Float.valueOf(lineStr[0])).floatValue();
			precursorIntensity = 0;
			precursorCharge = (Integer.valueOf(lineStr[1])).intValue();
			precursorMass = (float) calculateMassForCharge(0, 1, precursorMZ);
		} else {
		}
		tandemDict.put("PrecursorMass", precursorMass);
		tandemDict.put("PrecursorCharge", precursorCharge);
		tandemDict.put("PrecursorIntensity", precursorIntensity);
		while ((line = reader.readLine()) != null) {
			lineStr = line.split(" ");
			if (lineStr.length > 1) {
				massString = lineStr[0];
				intensityStr = lineStr[1];
				mOverZ = (Float.valueOf(massString)).floatValue();
				intensity = (Float.valueOf(intensityStr)).floatValue();
				if ((mOverZ > 55) && (mOverZ <= precursorMass)) {
					masses.add(mOverZ);
					intensities.add(intensity);
				}
			}
		}
		tandemDict.put("Masses", masses);
		tandemDict.put("Intensities", intensities);
		return tandemDict;

	}

	/********************************************************************
	 * If we use the top 100 high intense peaks the masses are no longer in
	 * order. Order them acccording to their masses (the way algorithm works)
	 ********************************************************************/
	public static void sortPeaksAsMasses(Hashtable tandemData) {
		Vector masses = (Vector) tandemData.get("Masses");
		Vector intensities = (Vector) tandemData.get("Intensities");
		int i, massCount = masses.size(); // NSNumber *mass;
		float currentMass, prevMass, currentIntensity, tempIntensity, tempMass;
		int Pos;
		float maxMass;
		for (i = 0; i < massCount; i++) {
			prevMass = (Float) masses.get(i);
			Pos = i;
			for (int j = i + 1; j < massCount; j++) {
				maxMass = (Float) masses.get(Pos);
				currentMass = (Float) masses.get(j);
				if (currentMass < maxMass) {
					Pos = j;

				}
			}
			tempIntensity = (Float) intensities.get(Pos);
			intensities.setElementAt(intensities.get(i), Pos);
			intensities.setElementAt(tempIntensity, i);
			tempMass = (Float) masses.get(Pos);
			masses.setElementAt(masses.get(i), Pos);
			masses.setElementAt(tempMass, i);

		}
	}

	public static void write(String s) {
	}

	public static void writeln(String s) {
	}

}
