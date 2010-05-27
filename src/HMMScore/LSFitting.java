package HMMScore;

import java.util.*;
import java.io.*;

public class LSFitting {


	public static void calculateFalsePositive(Hashtable statistics,
			Hashtable results) {
		Vector score = new Vector();
		Vector survival = new Vector();
		Vector histX = (Vector) statistics.get("Bins");

		Vector histY = (Vector) statistics.get("PeptideCount");
		int numOfPeptide = (Integer) statistics.get("TotalPeptideCount");
		score = histX;
		survival = calculateSurvival(histY);
		
		/*
		 * Find all of the values in the survival graph that 
		 * are greater than or equal to surMaxValue.
		 */
		int i = 0;
		int number1;
		double surMaxValue = 0.10;
		while ((Double) survival.get(i) > surMaxValue
				&& i < (survival.size() - 1))
			i++;
		number1 = i;
		
		/*
		 * find the first spot after the surMaxValue index where
		 * the histogram goes to 0
		 */
		int number2 = histY.size();
		for (i = number1; i < histY.size(); i++) {
			if (((Integer) histY.get(i)).intValue() == 0) {
				number2 = i;
				break;
			}
		}

		/*
		 * Now cull our lists based on these two indicies
		 */
		for (i = score.size() - 1; i >= number2; i--) {
			score.remove(i);
			survival.remove(i);
		}
		for (i = 0; i < number1; i++) {
			score.remove(0);
			survival.remove(0);
		}
		
		if (survival.size() > 5 && score.size() > 5) {
			fitData(score, survival, numOfPeptide, results);
		} else {
			results.put("FalsePositiveScore", new Vector());
		}

	}


	/************************************************************************************************/
	public static void fitData(Vector bins, Vector survival, int numPeptide,
			Hashtable results) {
		int i;
		double windowScore, falsePositive;
	
		Vector FPScore = new Vector();
		Vector scores = (Vector) results.get("Scores");
		int count = scores.size();
		survival = logSurvival(survival);
		double parameterM = parameterM(bins, survival);
		double parameterB = parameterB(bins, survival, parameterM);
		for (i = 0; i < count; i++) {
			double currentScore = (Double) scores.get(i);
	
			falsePositive = parameterM * currentScore + parameterB;
			falsePositive = Math.exp(falsePositive);
			falsePositive = numPeptide * falsePositive;
			FPScore.add(falsePositive);
		}
		results.put("FalsePositiveScore", FPScore);
	}


	public static Vector calculateSurvival(Vector histY) {
		int i, count = histY.size(), totalWindowCount = 0;
		int currentCount;
		double currentProbability, currentSurvival;
		Vector survival = new Vector();
		Vector prob = new Vector();
		for (i = 0; i < count; i++) {
			currentCount = (Integer) histY.get(i);
			totalWindowCount += currentCount;
		}
		for (i = 0; i < count; i++) {
			currentCount = (Integer) histY.get(i);
			currentProbability = (double) currentCount / totalWindowCount;
			prob.add(currentProbability);
		}
		for (i = 0; i < count; i++) {
			currentSurvival = calcSurvival(prob, i);
			survival.add(currentSurvival);
		}
		return survival;
	}

	public static double calcSurvival(Vector prob, int i) {
		int count = prob.size(), index;
		double survival = 0.0;
		for (index = i; index < count; index++)
			survival += (Double) prob.get(index);
		return survival;
	}


	public static Vector logSurvival(Vector score) {
		int i, count = score.size();
		double currentScore;
		double lScore;
		Vector logScore = new Vector();
		for (i = 0; i < count; i++) {
			currentScore = (Double) score.get(i);
			lScore = Math.log(currentScore);
			logScore.add(lScore);
		}
		return logScore;
	}

	/* This calculates b in y=mx+b equation */
	public static double parameterB(Vector xValues, Vector yValues,
			double parameterM) {
		double parameterB;
		double temp1, temp2;
		int count = xValues.size(), i;

		temp1 = temp2 = 0.0;
		for (i = 0; i < count; i++) {
			temp1 += (Integer) xValues.get(i);
		}

		for (i = 0; i < count; i++) {
			temp2 += (Double) yValues.get(i);
		}
		parameterB = (1.0 / count) * (temp2 - parameterM * temp1);
		return parameterB;
	}

	/******************************************************************************
	 *Now I will implement using the equation y=mx+b. In this case two
	 * parameters m and b are needed
	 ********************************************************************************/
	public static double parameterM(Vector xValues, Vector yValues) {
		double numerator1, numerator2, denomenator1, denomenator2;
		double numerator = 0.0, denomenator = 0.0;
		double temp1 = 0.0, temp2 = 0.0, temp = 0.0;
		double parameterM;
		int count = xValues.size(), i;
		for (i = 0; i < count; i++) {
			temp1 = (Integer) xValues.get(i);
			temp2 = (Double) yValues.get(i);
			temp += (temp1 * temp2);
		}
		numerator1 = count * (temp);
		temp1 = temp2 = 0.0;
		for (i = 0; i < count; i++) {
			temp1 += (Integer) xValues.get(i);
		}
		for (i = 0; i < count; i++) {
			temp2 += (Double) yValues.get(i);
		}
		numerator2 = temp1 * temp2;
		numerator = numerator1 - numerator2;
		temp1 = temp2 = 0.0;
		for (i = 0; i < count; i++) {
			temp1 = (Integer) xValues.get(i);
			temp2 += (temp1 * temp1);
		}
		denomenator1 = count * temp2;

		temp1 = 0.0;
		for (i = 0; i < count; i++) {
			temp1 += (Integer) xValues.get(i);
		}
		denomenator2 = (temp1 * temp1);
		denomenator = denomenator1 - denomenator2;
		parameterM = numerator / denomenator;
		return parameterM;
	}

}
