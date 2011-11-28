package Peppy;

import java.util.ArrayList;
import java.util.Collections;


public class ProteinDigestionServer {
	
	ArrayList<Protein> proteins;

	/* every protein produces a "cluster" of peptides */
	ArrayList<ArrayList<Peptide>> peptidesClusters;
	
	ArrayList<Thread> threads;
	
	/* this is how we keep track of which protein to give out next */
	private int proteinIndex = 0;
	
	/* why this simply isn't Properties.numberOfThreads is because we may have less proteins than that number */
	private int numberOfThreads;
	

	public ProteinDigestionServer(ArrayList<Protein> proteins) {
		this.proteins = proteins;	
		peptidesClusters = new ArrayList<ArrayList<Peptide>>(proteins.size());
		
		//here we make sure we don't use more threads than we have proteins
		numberOfThreads = Properties.numberOfThreads;
		if (Properties.numberOfThreads > proteins.size()) numberOfThreads = proteins.size();
		threads = new ArrayList<Thread>(numberOfThreads);
		
		//spawn new threads as needed
		for (int threadNumber = 0; threadNumber < numberOfThreads; threadNumber++) {
			ProteinDigestionThread digester = new ProteinDigestionThread(getNextProtein(), this);
			Thread thread = new Thread(digester);
			thread.start();
			threads.add(thread);	
		}
	}
	
	/**
	 * 
	 * @param peptidesFromOneProtein
	 * @return the next Protein on the list.  Null if there are no more.
	 */
	public synchronized Protein getNextProtein(ArrayList<Peptide> peptidesFromOneProtein) {
		peptidesClusters.add(peptidesFromOneProtein);
		return getNextProtein();
	}
	

	public synchronized Protein getNextProtein() {
		Protein out =  null;
		if (proteinIndex < proteins.size()) {
			out = proteins.get(proteinIndex);
			proteinIndex++;
		}
		return out;
	}
	
	/**
	 * First we wait for all of  our  threads to finish.
	 * @return accumulated 
	 */
	public ArrayList<Peptide> getPeptides() {
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
				e.printStackTrace();
			}
		}
		//calculate size of combined ArrayLists
		int size = 0;
		for (ArrayList<Peptide> peptideCluster: peptidesClusters) {
			size += peptideCluster.size();
		}
		//combine matches into result and return
		ArrayList<Peptide> out = new ArrayList<Peptide>(size);
		for (ArrayList<Peptide> peptideCluster: peptidesClusters) {
			out.addAll(peptideCluster);
		}
		/* sort the peptides */
		Collections.sort(out);
		return out;
	}

}
