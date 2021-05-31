package Experimental;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import Peppy.U;

public class SpectrumScoreStats {
	
	public static void main(String args[] ) {
		try {
			BufferedReader br = new BufferedReader(new FileReader("spectrum to score scatter.txt"));
			String line = br.readLine();
			
			int lineCount = 0;
			int greaterThanCount = 0;
			double maxScore = 0;
			double maxRatio = 0;
			while (line != null) {
				lineCount++;
				String [] chunks = line.split("\t");
				double peppyScore = Double.parseDouble(chunks[0]);
				double spectrumScore = Double.parseDouble(chunks[1]);
//				if (spectrumScore - peppyScore > 14.5) {
				if (peppyScore < 24 && peppyScore < spectrumScore) {
					greaterThanCount++;
					if (peppyScore > maxScore) maxScore = peppyScore;
				}
				
				/* ratio */
				double ratio = spectrumScore / peppyScore;
				if (ratio > maxRatio) maxRatio = ratio;
				line = br.readLine();
				
			}
			double percent = (double) greaterThanCount / lineCount;
			U.p(maxRatio);
			U.p(percent);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
