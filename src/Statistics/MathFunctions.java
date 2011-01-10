package Statistics;

import Utilities.U;


public class MathFunctions {
	
	private static final int length = 150;
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
				double total = 0;
				double probability;
				for (int i = k; i <= n; i++) {
					probability = MathFunctions.cachedNChooseK(n, i);
					probability *= Math.pow(p, i);
					probability *= Math.pow(1 - p, n - i);
					total += probability;
				}
				binomialProb50[n][k] = total;
			}
		}
		
		U.p("CachedMathFunctions initialized");
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
		double total = 0;
		double probability;
		for (int i = k; i <= n; i++) {
			probability = MathFunctions.cachedNChooseK(n, i);
			probability *= Math.pow(p, i);
			probability *= Math.pow(1 - p, n - i);
			total += probability;
		}
		return total;
	}
	
//	public static void main(String args[]) {
//		U.p("init");
//		double p = 0.5;
//
//		int n = 4;
//		int k = 1;
//		double total = 0;
//		double probability;
//		for (int i = k; i <= n; i++) {
//			probability = CachedMathFunctions.nChooseK(n, i);
//			probability *= Math.pow(p, i);
//			probability *= Math.pow(1 - p, n - i);
//			U.p(probability);
//			total += probability;
//		}
//		binomialProb50[n][k] = total;
//
//	}
	

}
