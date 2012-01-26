package Peppy;

import Utilities.U;


/**
 * Contains information about amino acids
 * Statically defined amino acids as bytes
 * amino acid masses
 * @author Brian Risk
 *
 */
public class AminoAcids {
	
	private static char [] acids = {'.', 'A', 'C',	'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
	
	private static byte tracker = 0;
	public static final byte STOP = tracker++;
	public static final byte A = tracker++;
	public static final byte C = tracker++;
	public static final byte D = tracker++;
	public static final byte E = tracker++;
	public static final byte F = tracker++;
	public static final byte G = tracker++;
	public static final byte H = tracker++;
	public static final byte I = tracker++;
	public static final byte K = tracker++;
	public static final byte L = tracker++;
	public static final byte M = tracker++;
	public static final byte N = tracker++;
	public static final byte P = tracker++;
	public static final byte Q = tracker++;
	public static final byte R = tracker++;
	public static final byte S = tracker++;
	public static final byte T = tracker++;
	public static final byte V = tracker++;
	public static final byte W = tracker++;
	public static final byte Y = tracker++;
	
	
	private static double [] weightsMono = new double[acids.length];
	private static double [] weightsAverage = new double[acids.length];
	
	static {
		init();
	}
	
	public static void init(){
		weightsMono[STOP] = 0.0;
		weightsMono[A] = 71.03711;
		weightsMono[C] = 103.00919;
		weightsMono[D] = 115.02694;
		weightsMono[E] = 129.04259;
		weightsMono[F] = 147.06841;
		weightsMono[G] = 57.02146;
		weightsMono[H] = 137.05891;
		weightsMono[I] = 113.08406;
		weightsMono[K] = 128.09496;
		weightsMono[L] = 113.08406;
		weightsMono[M] = 131.04049;
		weightsMono[N] = 114.04293;
		weightsMono[P] = 97.05276;
		weightsMono[Q] = 128.05858;
		weightsMono[R] = 156.10111;
		weightsMono[S] = 87.03203;
		weightsMono[T] = 101.04768;
		weightsMono[V] = 99.06841;
		weightsMono[W] = 186.07931;
		weightsMono[Y] = 163.06333;
		
		weightsAverage[STOP] = 0.0;
		weightsAverage[A] = 71.0788;
		weightsAverage[C] = 103.1448;
		weightsAverage[D] = 115.0886;
		weightsAverage[E] = 129.1155;
		weightsAverage[F] = 147.1766;
		weightsAverage[G] = 57.052;
		weightsAverage[H] = 137.1412;
		weightsAverage[I] = 113.1595;
		weightsAverage[K] = 128.1742;
		weightsAverage[L] = 113.1595;
		weightsAverage[M] = 131.1986;
		weightsAverage[N] = 114.1039;
		weightsAverage[P] = 97.1167;
		weightsAverage[Q] = 128.1308;
		weightsAverage[R] = 156.1876;
		weightsAverage[S] = 87.0782;
		weightsAverage[T] = 101.1051;
		weightsAverage[V] = 99.1326;
		weightsAverage[W] = 186.2133;
		weightsAverage[Y] = 163.176;
		
		if(Properties.useIsotopeLabeling) {
			U.p("amino acid weights adjusted for labeling");
			weightsMono[K] = 136.109159;
			weightsAverage[K] = 136.1169;
			weightsMono[R] = 166.109379;
			weightsAverage[R] = 166.1172;	
		}
	}
	
	
	public static double getWeightMono(byte acid) {
		return weightsMono[acid];
	}
	
	public static double getWeightMono(char acid) {
		return getWeightMono(getAminoAcidByte(acid));
	}

	
	public static double getWeightAverage(byte acid) {
		return weightsAverage[acid];
	}
	
	public static double getWeightAverage(char acid) {
		return getWeightAverage(getAminoAcidByte(acid));
	}
	
	
	public static char getAminoAcidChar(byte acid) {
		if (acid == -1) return 'X';
		return acids[acid];
	}

	
	public static byte getAminoAcidByte(char acid) {
		if (acid == '.') return STOP;
		if (acid == 'A') return A;
		if (acid == 'C') return C;
		if (acid == 'D') return D;
		if (acid == 'E') return E;
		if (acid == 'F') return F;
		if (acid == 'G') return G;
		if (acid == 'H') return H;
		if (acid == 'I') return I;
		if (acid == 'K') return K;
		if (acid == 'L') return L;
		if (acid == 'M') return M;
		if (acid == 'N') return N;
		if (acid == 'P') return P;
		if (acid == 'Q') return Q;
		if (acid == 'R') return R;
		if (acid == 'S') return S;
		if (acid == 'T') return T;
		if (acid == 'V') return V;
		if (acid == 'W') return W;
		if (acid == 'Y') return Y;
		return STOP;
//		throw new Error("Invalid amino acid character: " + acid);
	}
	
	public static boolean isValid(byte acid) {
		if (acid == STOP) return true;
		if (acid == A) return true;
		if (acid == C) return true;
		if (acid == D) return true;
		if (acid == E) return true;
		if (acid == F) return true;
		if (acid == G) return true;
		if (acid == H) return true;
		if (acid == I) return true;
		if (acid == K) return true;
		if (acid == L) return true;
		if (acid == M) return true;
		if (acid == N) return true;
		if (acid == P) return true;
		if (acid == Q) return true;
		if (acid == R) return true;
		if (acid == S) return true;
		if (acid == T) return true;
		if (acid == V) return true;
		if (acid == W) return true;
		if (acid == Y) return true;
		return false;
	}
	
	
	public static byte [] getByteArrayForString(String acidString) {
		byte [] out = new byte[acidString.length()];
		for (int i = 0; i < acidString.length(); i++) {
			out[i] = getAminoAcidByte(acidString.charAt(i));
		}
		return out;
	}
	
	
	public static String getStringForByteArray(byte [] acidArray) {
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < acidArray.length; i++) {
			sb.append(getAminoAcidChar(acidArray[i]));
		}
		return sb.toString();
	}
	
	public static int getNumberOfAminoAcids() {
		return acids.length;
	}
	
}
