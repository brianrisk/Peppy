package Math;

import java.util.ArrayList;




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
		return Math.exp(n * Math.log(n) - n + 1);
	}
	
	public static double getCachedBinomialProbability50(int n, int k) {
		return binomialProb50[n][k];
	}

	public static double getBinomialProbability(int n, int k, double p) {
		if (k <= n * p) return 1;
		double total = 0.0;
		double probability;
		for (int i = k; i <= n; i++) {
			probability = MathFunctions.cachedNChooseK(n, i);
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
