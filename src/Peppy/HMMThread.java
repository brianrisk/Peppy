package Peppy;

import java.util.ArrayList;
import java.util.Collections;

public class HMMThread implements Runnable {
	
	Match match;
	HMMScorer scorer;
	
	public HMMThread(Match match, HMMScorer scorer) {
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
