package Peppy;

import java.util.*;

public class HMMScorer {

	private ArrayList<Match> matches;
	private int numberOfThreads = Properties.numberOfThreads;
	private int matchIndex = 0;
	private ArrayList<Thread> threads = new ArrayList<Thread>(Properties.numberOfThreads);
	
	public HMMScorer(ArrayList<Match> matches) {
		this.matches = matches;
		
		//here we make sure we don't use more threads than we have matches
		if (numberOfThreads > matches.size()) numberOfThreads = matches.size();
	}
	
	public void score() {
		//spawn new threads as needed
		for (int threadNumber = 0; threadNumber < numberOfThreads; threadNumber++) {
			HMMThread scorer = new HMMThread(getNextMatch(), this);
			Thread thread = new Thread(scorer);
			thread.start();
			threads.add(thread);	
		}
		//Now wait for them all to finish
		for (Thread thread: threads) {
			try {
				thread.join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
	}

	public synchronized Match getNextMatch() {
		Match out =  null;
		if (matchIndex < matches.size()) {
			out = matches.get(matchIndex);
			matchIndex++;
		}
		return out;
	}
	
	

}
