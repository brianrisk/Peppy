package HMMScore;

public class Defines {

	public final static int numberOfIons = 12;
	public final static int bin_count = 10;
	public final static int numberOfAminoAcid = 20;

	public final static int B_ION = 0;
	public final static int Y_ION = 1;
	public final static int A_ION = 2;
	public final static int J_ION = 3;
	public final static int B17_ION = 4;
	public final static int Y17_ION = 5;
	public final static int B18_ION = 6;
	public final static int Y18_ION = 7;
	public final static int IMM_ION = 8;
	public final static int INTERNAL_ION = 9;
	public final static int INTERNAL17 = 10;
	public final static int INTERNAL18 = 11;

	/*
	 * public final static int C_ION =13; public final static int Z_ION =14;
	 * public final static int X_ION =15; public final static int A18_ION= 16;
	 * public final static int A17_ION = 17;
	 */

	public final static double OXYGEN_MONO = 15.99491463;
	public final static double OXYGEN_AVE = 15.99940494;
	public final static double NITROGEN_MONO = 14.00307400;
	public final static double NITROGEN_AVE = 14.00674309;
	public final static double HYDROGEN_MONO = 1.00782504;
	public final static double HYDROGEN_AVE = 1.00794076;
	public final static double CARBON_MONO = 12.00000000;
	public final static double CARBON_AVE = 12.01073590;
	public final static double PHOSPHORUS_MONO = 30.973762;
	public final static double SULPHUR_MONO = 31.972070698;
	public final static double AMINO_ACID_NUMBER = 20;
	public final static double WATER_MONO = (2 * HYDROGEN_MONO) + OXYGEN_MONO;
	public final static double AMONIA_MONO = (NITROGEN_MONO + 3 * HYDROGEN_MONO);
	public final static double WATER_AVERAGE = OXYGEN_AVE + 2 * HYDROGEN_AVE;
	public final static double AMONIA_AVERAGE = NITROGEN_AVE + 3 * HYDROGEN_AVE;
	public final static double CO = CARBON_MONO + OXYGEN_MONO; // mass of carbon
																// monoxide
	public final static int MONOISOTOPIC = 0;
	public final static int AVERAGE = 1;
	public final static double THRESHOLD = 0.5;
	public final static int TOTAL_PEAK_COUNT = 100;
	public final static double PEP_TOLERANCE = 2.0;

	public final static int BIN1 = 0;
	public final static int BIN2 = 1;
	public final static int BIN3 = 2;
	public final static int BIN4 = 3;
	public final static int BIN5 = 4;
	public final static int BIN6 = 5;
	public final static int BIN7 = 6;
	public final static int BIN8 = 7;
	public final static int BIN9 = 8;
	public final static int BIN10 = 9;

	public static final int A = 0;
	public static final int R = 1;
	public static final int N = 2;
	public static final int D = 3;
	public static final int C = 4;
	public static final int E = 5;
	public static final int Q = 6;
	public static final int G = 7;
	public static final int H = 8;
	public static final int I = 9;
	public static final int L = 10;
	public static final int K = 11;
	public static final int M = 12;
	public static final int F = 13;
	public static final int P = 14;
	public static final int S = 15;
	public static final int T = 16;
	public static final int W = 17;
	public static final int Y = 18;
	public static final int V = 19;
}
