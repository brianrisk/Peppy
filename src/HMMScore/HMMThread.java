package HMMScore;

import Peppy.Match;


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
