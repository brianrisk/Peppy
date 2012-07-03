package Peppy;

import java.util.ArrayList;

/**
 * Copyright 2012, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class VariableModificationServer {
	
	ArrayList<Peptide> peptides;
	
	ArrayList<Thread> threads;
	
	//this is how we keep track of which Spectrum to give out next
	private int spectrumMatchesIndex = 0;
	
	/* why this simply isn't Properties.numberOfThreads is because we may have less spectra than that number */
	private int numberOfThreads;
	
//	/**
//	 * 
//	 * @param peptides
//	 * @param spectra
//	 * @param matches the ArrayList where we store the best matches
//	 */
//	public ScoringServer(ArrayList<Peptide> peptides, ArrayList<MatchesSpectrum> spectraMatches) {
//		this.spectraMatches = spectraMatches;	
//		
//		//here we make sure we don't use more threads than we have spectra
//		numberOfThreads = Properties.numberOfThreads;
//		if (numberOfThreads > spectraMatches.size()) numberOfThreads = spectraMatches.size();
//		threads = new ArrayList<Thread>(numberOfThreads);
//		
//		//spawn new threads as needed
//		for (int threadNumber = 0; threadNumber < numberOfThreads; threadNumber++) {
//			ScoringThread scorer = new ScoringThread(getNextSpectrumMatches(), peptides, this);
//			Thread thread = new Thread(scorer);
//			thread.start();
//			threads.add(thread);	
//		}
//	}
	
	/**
	 * takes no parameters because it is used at the start when there are no matches to incorporate.
	 * @return
	 */
	public synchronized MatchesSpectrum getNextPeptide() {
		MatchesSpectrum out =  null;
//		if (spectrumMatchesIndex < spectraMatches.size()) {
//			out = spectraMatches.get(spectrumMatchesIndex);
//			spectrumMatchesIndex++;
//		}
		return out;
	}
	
	/**
	 * First we wait for all of  our threads to finish.
	 * "Yo threads, I'm happy for you and I'ma let you finish"
	 * @return accumulated 
	 */
	public void findMatches() {
		boolean going = true;
		while (going) {
			for (Thread thread: threads) {
				going = thread.isAlive();
				//if at least one thread is going, break out of this for loop
				if (going) break;
			}
			//Sleep for a bit and wait for the threads to finish.
			try {
				Thread.sleep(500);
			} catch (InterruptedException e) {
				U.p("getMatches in ScoringEngine thread was interrupted!");
				e.printStackTrace();
			}
		}
		
	}

}
