package Math;

import java.util.ArrayList;

import Utilities.U;




public class MathFunctions {
	
	private static final int length = 200;
	private static double [] logs = new double[length];
	private static double [] factorial = new double[length];
	private static double [][] nChooseK = new double [length][length];
	private static double [][] binomialProb50 = new double [length][length];
	
	
	static {
		//logs
		for (int i = 0; i < length; i++) {
			logs[i] = Math.log(i);
		}
		
		//factorial
		double nFactorial = 1;
		factorial[0] = 1;
		for (int n = 1; n < length; n++) {
			nFactorial *= n;
			factorial[n] = nFactorial;
		}
		
		//n choose k
		for (int n = 1; n < length; n++) {
			for (int k = 0; k <= n; k++) {
				nChooseK[n][k] = factorial[n] / (factorial[k] * factorial[n -k]);
			}
		}
		
		//binomial prob 50
		double p = 0.5;
		for (int n = 1; n < length; n++) {
			for (int k = 0; k <= n; k++) {
				binomialProb50[n][k] = getBinomialProbability(n, k, p);
			}
		}
		
	}
	
	public static void main(String args[]) {
//		int max = 20;
//		for (int i = 1; i < max; i++) {
//			U.p(approximateBinomialProbability(max, i, 0.5) / getCachedBinomialProbability50(max, i));
////			U.p(i + ": " + approximateBinomialProbability(max, i, 0.5));
////			U.p(cachedNChooseK(max, i) / approximateNChooseK(max, i));
////			U.p(i + ": " + approximateNChooseK(max, i) / cachedNChooseK(max, i));
//		}
//		U.p(approximateBinomialProbability(527, 10, 0.014084507042253521));
	}
	
	public static double cachedLog(int n) {
		return logs[n];
	}
	
	public static double cachedNChooseK(int n, int k) {
		return nChooseK[n][k];
	}
	
	public static double approximateNChooseK(int n, int k) {
		return approximateFactorial(n) / (approximateFactorial(k) * approximateFactorial(n - k));
	}
	
	public static double approximateFactorial(int n) {
		if (n == 0) return 1;
		double out = Math.sqrt(2.0 * Math.PI * n);
		out *= Math.pow((double) n / Math.E, n);
		return out;
	}
	
	public static double approximateBinomialProbability(int n, int k, double p) {
		if (k <= n * p) return 1;
		double total = 0.0;
		double mean = n * p;
		double variance = mean * ( 1.0 - p);
		for (int i = k; i <= n; i++) {
			total += gaussian(i, mean, variance);
		}
		return total;
	}
	
	public static double approximateSingleBinomialProbability(int n, int k, double p) {
		double mu = n * p;
		double sigma = Math.sqrt(mu * ( 1.0 - p));
		return (k - mu) / sigma;
	}
	
	public static double gaussian(double x, double mean, double variance) {
		double exponent = - (x - mean) * (x - mean);
		exponent /= 2.0 * variance;
		double out = 1.0 / (Math.sqrt(2 * Math.PI * variance));
		out *= Math.exp(exponent);
		return out;
		
	}
	
	public static double getCachedBinomialProbability50(int n, int k) {
		return binomialProb50[n][k];
	}

	public static double getBinomialProbability(int n, int k, double p) {
		if (k <= n * p) return 1;
		double total = 0.0;
		double probability;
		for (int i = k; i <= n; i++) {
			probability = cachedNChooseK(n, i);
			probability *= Math.pow(p, i);
			probability *= Math.pow(1 - p, n - i);
			total += probability;
		}
		return total;
	}
	
	/**
	 * Boolean search to locate the first peptide in the SORTED list of peptides that has
	 * a mass greater than the "mass" parameter.
	 * 
	 * CAUTION:  this method is not perfect due to rounding error.  However, returns
	 * very good ballpark.
	 * 
	 * @param values
	 * @param value
	 * @return
	 */
	public static int findFirstIndexGreater(ArrayList<? extends HasValue> values, double value) {
		int index = values.size() / 2;
		int increment = index / 2;
		while (increment > 0) {
			if (values.get(index).getValue() > value) {index -= increment;}
			else {index += increment;}
			increment /= 2;
		}
		return index;
	}
	

}
