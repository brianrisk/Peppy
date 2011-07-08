package Peppy;


public class ProteinDigestionThread implements Runnable {
	
	Protein protein;
	ProteinDigestionServer proteinDigestionServer;
	
	/**
	 * @param proteins
	 * @param spectrum
	 */
	public ProteinDigestionThread(Protein protein, ProteinDigestionServer proteinDigestionServer) {
		this.protein = protein;
		this.proteinDigestionServer = proteinDigestionServer;
	}
	

	public void run() {
		while (protein != null) {
			//return results, get new task
			protein = proteinDigestionServer.getNextProtein(protein.getPeptides());
		}
	}

	
	

}
