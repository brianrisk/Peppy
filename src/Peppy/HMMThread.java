package Peppy;

import java.util.ArrayList;
import java.util.Collections;

public class HMMThread implements Runnable {
	
	SpectrumPeptideMatch match;
	HMMScorer scorer;
	
	public HMMThread(SpectrumPeptideMatch match, HMMScorer scorer) {
		this.match = match;
		this.scorer = scorer;
	}

	public void run() {
		while (match != null) {
			match.calculateHMM();
			match = scorer.getNextMatch();
		}
	}

}
