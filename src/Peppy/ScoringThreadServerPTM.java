package Peppy;
import java.util.ArrayList;
import java.util.Collections;

import Utilities.U;


public class ScoringThreadServerPTM {
	
	ArrayList<Peptide> peptides;
	/*
	 * Having a ArrayList of ArrayLists may look overly complicated, but I am doing this for a reason.
	 * The easier way is to just have a ArrayList of SpectrumPeptideMatch objects and each time 
	 * getNextSpectrum we take the result we're given and "addAll".  The unfortunate thing is
	 * that since getNextSpectrum is synchronized it makes all the other possibly dozens of
	 * threads wait.  Let's not make them wait.  Let's let them do their work as fast as they
	 * can and then merge all the match ArrayLists together at the end.
	 */
	ArrayList<ArrayList<Match_IMP_VariMod>> matches;
	
	ArrayList<Thread> threads = new ArrayList<Thread>(Properties.numberOfThreads);
	
	//this is how we keep track of which Spectrum to give out next
	private int peptideIndex = 0;
	private int numberOfThreads;
	
	/**
	 * 
	 * @param peptides
	 * @param spectra
	 * @param matches the ArrayList where we store the best matches
	 */
	public ScoringThreadServerPTM(ArrayList<Peptide> peptides, ArrayList<Spectrum> spectra) {
		this.peptides = peptides;
		matches = new ArrayList<ArrayList<Match_IMP_VariMod>>(peptides.size());
		
		//making sure lists are sorted
		Collections.sort(spectra);
		Collections.sort(peptides);
		
		//here we make sure we don't use more threads than we have spectra
		numberOfThreads = Properties.numberOfThreads;
		if (numberOfThreads > peptides.size()) numberOfThreads = peptides.size();
		
		//spawn new threads as needed
		for (int threadNumber = 0; threadNumber < numberOfThreads; threadNumber++) {
			ScoringThreadPTM scorer = new ScoringThreadPTM(getNextPeptide(), spectra, this);
			Thread thread = new Thread(scorer);
			thread.start();
			threads.add(thread);	
		}
	}
	
	/**
	 * 
	 * @param matchesForOneSpectrum
	 * @return the next Spectrum on the list.  Null if there are no more.
	 */
	public synchronized Peptide getNextPeptide(ArrayList<Match_IMP_VariMod> matchesForOneSpectrum) {
		matches.add(matchesForOneSpectrum);
		return getNextPeptide();
	}
	
	/**
	 * takes no parameters because it is used at the start when there are no matches to incorporate.
	 * @return
	 */
	public synchronized Peptide getNextPeptide() {
		Peptide out =  null;
		if (peptideIndex < peptides.size()) {
			out = peptides.get(peptideIndex);
			peptideIndex++;
		}
		return out;
	}
	
	/**
	 * First we wait for all of  our scoring threads to finish.
	 * "Yo ScoringThreads, I'm happy for you and I'ma let you finish"
	 * @return accumulated 
	 */
	public ArrayList<Match_IMP_VariMod> getMatches() {
		boolean going = true;
		while (going) {
			for (int threadNumber = 0; threadNumber < numberOfThreads; threadNumber++) {
				Thread thread = threads.get(threadNumber);
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
		//calculate size of combined ArrayLists
		int size = 0;
		for (ArrayList<Match_IMP_VariMod> matchCluster: matches) {
			size += matchCluster.size();
		}
		//combine matches into result and return
		ArrayList<Match_IMP_VariMod> out = new ArrayList<Match_IMP_VariMod>(size);
		for (ArrayList<Match_IMP_VariMod> matchCluster: matches) {
			out.addAll(matchCluster);
		}
		return out;
	}

}
