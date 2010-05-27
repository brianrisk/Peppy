package HMMScore;

import java.io.*;
import java.util.*;

public class HMM_SCORE {
	public static double[] AAMasses = new double[200];

	public static void setAAMonoMasses() {
		AAMasses['A'] = 71.03711;
		AAMasses['C'] = 103.00919;
		AAMasses['D'] = 115.02694;
		AAMasses['L'] = 113.08406;
		AAMasses['I'] = 113.08406;
		AAMasses['K'] = 128.09496;
		AAMasses['N'] = 114.04293;
		AAMasses['T'] = 101.04768;
		AAMasses['R'] = 156.10111;
		AAMasses['Y'] = 163.06333;
		AAMasses['V'] = 99.06841;
		AAMasses['Q'] = 128.05858;
		AAMasses['E'] = 129.04259;
		AAMasses['G'] = 57.02146;
		AAMasses['H'] = 137.05891;
		AAMasses['M'] = 131.04049;
		AAMasses['F'] = 147.06841;
		AAMasses['P'] = 97.05276;
		AAMasses['S'] = 87.03203;
		AAMasses['W'] = 186.07931;
	}
}
