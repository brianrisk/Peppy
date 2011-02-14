package HMMScore;
import java.util.ArrayList;
import java.util.Vector;

import Peppy.Peak;

public class EmitionProb {

	static float[][] massEmitMatrix = new float[Defines.numberOfIons][Defines.bin_count];
	static float[][] intensityEmitMatrix = new float[Defines.numberOfIons][Defines.bin_count];

	/*********************************************************************************/
	public static int getBinForMass(double mass, double PMass) {
		double ratio = mass / PMass;
		if (ratio <= 0.1) {
			return Defines.BIN1;
		}
		if (ratio > 0.1 & ratio <= 0.2) {
			return Defines.BIN2;
		}
		if (ratio > 0.2 & ratio <= 0.3) {
			return Defines.BIN3;
		}

		if (ratio > 0.3 & ratio <= 0.4) {
			return Defines.BIN4;
		}

		if (ratio > 0.4 & ratio <= 0.5) {
			return Defines.BIN5;
		}
		if (ratio > 0.5 & ratio <= 0.6) {
			return Defines.BIN6;
		}
		if (ratio > 0.6 & ratio <= 0.7) {
			return Defines.BIN7;
		}
		if (ratio > 0.7 & ratio <= 0.8) {
			return Defines.BIN8;
		}
		if (ratio > 0.8 & ratio <= 0.9) {
			return Defines.BIN9;
		}
		if (ratio > 0.9 & ratio <= 1.0) {
			return Defines.BIN10;
		}
		return Defines.BIN10;

	}

	/***************************************************************************/
	public static float[][] getIntensityEmitProbMat() {
		return intensityEmitMatrix;
	}

	/************************************************************************************/
	public static int getBinForIntensity(double intensity,
			Vector<Double> intensities) {
		int il, currentBin = 0;
		double intensityToCheck;
		while (currentBin < Defines.bin_count) {
			Vector<Double> currentBinIntensities = getBinIntensities(
					intensities, currentBin);
			for (il = 0; il < currentBinIntensities.size(); il++) {
				intensityToCheck = currentBinIntensities.get(il);
				if (Math.abs(intensity - intensityToCheck) < 0.01)
					return currentBin;
			}
			currentBin++;
		}
		return currentBin;
	}

	/*******************************************************************************************/
	public static Vector<Double> getBinIntensities(Vector<Double> intensities,
			int currentBin) {
		Vector<Double> retArray = new Vector<Double>();
		int ik = 0, totalIntensityCount = intensities.size();
		int currentBinCount = (int) (0.1 * totalIntensityCount + 0.2);// 10% in
																		// each
																		// bin
		if (currentBin == 0) {
			ik = 0;
			while (ik < currentBinCount) {
				retArray.add(intensities.get(ik));
				ik++;
			}
			return retArray;
		} else if (currentBin == 1) {
			ik = currentBinCount;
			while (ik < 2 * currentBinCount) {
				retArray.add(intensities.get(ik));
				ik++;
			}
			return retArray;
		}

		else if (currentBin == 2) {
			ik = 2 * currentBinCount;
			while (ik < 3 * currentBinCount) {
				retArray.add(intensities.get(ik));
				ik++;
			}
			return retArray;
		} else if (currentBin == 3) {
			ik = 3 * currentBinCount;
			while (ik < 4 * currentBinCount) {
				retArray.add(intensities.get(ik));
				ik++;
			}
			return retArray;
		} else if (currentBin == 4) {
			ik = 4 * currentBinCount;
			while (ik < 5 * currentBinCount) {
				retArray.add(intensities.get(ik));
				ik++;
			}
			return retArray;
		} else if (currentBin == 5) {
			ik = 5 * currentBinCount;
			while (ik < 6 * currentBinCount) {
				retArray.add(intensities.get(ik));
				ik++;
			}
			return retArray;
		} else if (currentBin == 6) {
			ik = 6 * currentBinCount;
			while (ik < 7 * currentBinCount) {
				retArray.add(intensities.get(ik));
				ik++;
			}
			return retArray;
		} else if (currentBin == 7) {
			ik = 7 * currentBinCount;
			while (ik < 8 * currentBinCount) {
				retArray.add(intensities.get(ik));
				ik++;
			}
			return retArray;
		} else if (currentBin == 8) {
			ik = 8 * currentBinCount;
			while (ik < 9 * currentBinCount) {
				retArray.add(intensities.get(ik));
				ik++;
			}
			return retArray;
		} else {
			ik = 9 * currentBinCount;
			while (ik < totalIntensityCount) {
				retArray.add(intensities.get(ik));
				ik++;
			}
			return retArray;
		}

	}

	/***********************************************************************************
	 * Delete some peaks if masses are less or greater than some certain values:
	 * Jainab
	 ***********************************************************************************/
	public static ArrayList<Peak> cleanPeakForHighMass(ArrayList<Peak> peaks,
			double precursorMass) {
		ArrayList<Peak> cleanPeaks = new ArrayList<Peak>();
		double cmOverZ;
		for (Peak peak : peaks) {
			cmOverZ = peak.getMass();
			if ((cmOverZ > 70) && (cmOverZ <= precursorMass)) {
				cleanPeaks.add(peak);
			}
		}
		return cleanPeaks;

	}

	/******************************************************************************************
	 * public static int getIonTypeForMass(double mass, Hashtable allTypeIons) {
	 * int i; Vector currentTypeIons;
	 * 
	 * currentTypeIons = (Vector) allTypeIons.get("Y_IONS"); //[allTypeIons
	 * objectForKey:@"Y_IONS"]; for (i=0; i< currentTypeIons.size() ; i++) { if
	 * (Math.abs( mass -(Double) currentTypeIons.get(i) ) < Properties.peakDifferenceThreshold )
	 * { return Defines.Y_ION; } }
	 * 
	 * currentTypeIons = (Vector)allTypeIons.get("B_IONS"); for (i=0; i<
	 * currentTypeIons.size() ; i++) { if (Math.abs( mass -(Double)
	 * currentTypeIons.get(i) ) < Properties.peakDifferenceThreshold ) { return Defines.B_ION; }
	 * } currentTypeIons =(Vector) allTypeIons.get("Y17_IONS"); for (i=0; i<
	 * currentTypeIons.size() ; i++) { if (Math.abs( mass -(Double)
	 * currentTypeIons.get(i) ) < Properties.peakDifferenceThreshold ) { return Defines.Y17_ION;
	 * } } currentTypeIons =(Vector) allTypeIons.get("Y18_IONS"); for (i=0; i<
	 * currentTypeIons.size() ; i++) { if (Math.abs( mass -(Double)
	 * currentTypeIons.get(i) ) < Properties.peakDifferenceThreshold ) { return Defines.Y18_ION;
	 * } } currentTypeIons =(Vector) allTypeIons.get("B17_IONS"); for (i=0; i<
	 * currentTypeIons.size() ; i++) { if (Math.abs( mass -(Double)
	 * currentTypeIons.get(i) ) < Properties.peakDifferenceThreshold ) { return Defines.B17_ION;
	 * } } currentTypeIons =(Vector) allTypeIons.get("B18_IONS"); for (i=0; i<
	 * currentTypeIons.size() ; i++) { if (Math.abs( mass -(Double)
	 * currentTypeIons.get(i) ) < Properties.peakDifferenceThreshold ) { return Defines.B18_ION;
	 * } } currentTypeIons =(Vector) allTypeIons.get("A_IONS"); for (i=0; i<
	 * currentTypeIons.size() ; i++) { if (Math.abs( mass -(Double)
	 * currentTypeIons.get(i) ) < Properties.peakDifferenceThreshold ) { return Defines.A_ION; }
	 * }
	 * 
	 * currentTypeIons =(Vector) allTypeIons.get("IMM_IONS"); for (i=0; i<
	 * currentTypeIons.size() ; i++) { if (Math.abs( mass -(Double)
	 * currentTypeIons.get(i) ) < Properties.peakDifferenceThreshold ) { return Defines.IMM_ION;
	 * } }
	 * 
	 * 
	 * currentTypeIons =(Vector) allTypeIons.get("INTERNAL"); for (i=0; i<
	 * currentTypeIons.size() ; i++) { if (Math.abs( mass -(Double)
	 * currentTypeIons.get(i) ) < Properties.peakDifferenceThreshold ) { return
	 * Defines.INTERNAL_ION; } } currentTypeIons =(Vector)
	 * allTypeIons.get("INTERNAL17"); for (i=0; i< currentTypeIons.size() ; i++)
	 * { if (Math.abs( mass -(Double) currentTypeIons.get(i) ) <
	 * Properties.peakDifferenceThreshold ) { return Defines.INTERNAL17; } } currentTypeIons
	 * =(Vector) allTypeIons.get("INTERNAL18"); for (i=0; i<
	 * currentTypeIons.size() ; i++) { if (Math.abs( mass -(Double)
	 * currentTypeIons.get(i) ) < Properties.peakDifferenceThreshold ) { return
	 * Defines.INTERNAL18; } }
	 */

	/***********************************************************************************
	 * Get the top 100 most intense peaks: Jainab
	 ***********************************************************************************/
	public static ArrayList<Peak> getHighIntensityPeaks(ArrayList<Peak> peaks) {
		int peakCount = peaks.size();
		int i, j;
		Peak tempPeak, currentPeak;
		ArrayList<Peak> retVector = new ArrayList<Peak>();

		for (i = 1; i < peakCount; i++) {
			currentPeak = peaks.get(i);
			j = i;
			while (j > 0
					&& peaks.get(j - 1).getIntensity() < currentPeak
							.getIntensity()) {
				tempPeak = peaks.get(j);
				peaks.set(i, peaks.get(j - 1));
				peaks.set(j - 1, tempPeak);
				j--;
			}
			peaks.set(j, currentPeak);
		}

		/*
		 * cout << " After sorting print the peaks  \n"; for (int j=0;
		 * j<peaks.size(); j++) { cout << "Print peaks here " << peaks[j].m_fM
		 * <<"\t" <<peaks[j].m_fI<<"\n"; }
		 */

		for (int ii = 0; ii < Defines.TOTAL_PEAK_COUNT; ii++) {
			retVector.add(peaks.get(ii));
		}

		return retVector;

	}

	/***************************************************************************/
	public static float[][] getMassEmitProbMat() {
		return massEmitMatrix;
	}

	/*************************************************************************************************************/
	public static Vector<Double> sortIntensities(ArrayList<Peak> peaks) {
		Vector<Double> intensities = new Vector<Double>();
		for (Peak peak : peaks) {
			intensities.add((double) peak.getIntensity());
		}
		int i, j, intensityCount = intensities.size();
		double tempIntensity, prevIntensity, maxIntensity, currentIntensity;
		int maxIntensityPos;
		for (i = 0; i < intensityCount - 1; i++) {
			prevIntensity = intensities.get(i);
			maxIntensityPos = i;
			for (j = i + 1; j < intensityCount; j++) {
				maxIntensity = intensities.get(maxIntensityPos);
				currentIntensity = intensities.get(j);
				if (currentIntensity > maxIntensity) {
					maxIntensityPos = j;

				}
			}
			tempIntensity = intensities.get(maxIntensityPos); // [
																// [mutableIntensities
																// objectAtIndex:maxIntensityPos]
																// floatValue];
			intensities.setElementAt(intensities.get(i), maxIntensityPos);
			intensities.setElementAt(tempIntensity, i);

		}
		return intensities;
	}

	public static void write(String s) {
	}

	public static void writeln(String s) {
	}

}
