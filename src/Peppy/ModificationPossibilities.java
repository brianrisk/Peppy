package Peppy;

import java.util.ArrayList;

import Utilities.U;

public class ModificationPossibilities {
	@SuppressWarnings("unchecked")
	ArrayList[] modificationsArray;
	
	public ModificationPossibilities() {
		modificationsArray = new ArrayList[AminoAcids.getNumberOfAminoAcids()];
		for (int i = 0; i < AminoAcids.getNumberOfAminoAcids(); i++) {
			modificationsArray[i] = new ArrayList<Modification>();
		}
		
		//TODO remove later
		//specific to Yanbao

		//Acetylation	K, R
		addModificationToAminoAcid("Acetylation", AminoAcids.K);
		addModificationToAminoAcid("Acetylation", AminoAcids.R);
		//ADP ribosylation	E
		addModificationToAminoAcid("ADP ribosylation", AminoAcids.E);
		//Biotinylation	K
		addModificationToAminoAcid("Biotinylation", AminoAcids.K);
		//Butyrylation	K
		addModificationToAminoAcid("Butyrylation", AminoAcids.K);
		//Dimethylation	K, R
		addModificationToAminoAcid("Dimethylation", AminoAcids.K);
		addModificationToAminoAcid("Dimethylation", AminoAcids.R);
		//Methylation	K, R
		addModificationToAminoAcid("Methylation", AminoAcids.K);
		addModificationToAminoAcid("Methylation", AminoAcids.R);
		//Oxidation	M
		addModificationToAminoAcid("Oxidation", AminoAcids.M);
		//Palmitoylation	K
		addModificationToAminoAcid("Palmitoylation", AminoAcids.K);
		//Phosphorylation	S, T, Y
		addModificationToAminoAcid("Phosphorylation", AminoAcids.S);
		addModificationToAminoAcid("Phosphorylation", AminoAcids.T);
		addModificationToAminoAcid("Phosphorylation", AminoAcids.Y);
		//Propionylation	K
		addModificationToAminoAcid("Propionylation", AminoAcids.K);
		//Sumoylation1	K
		addModificationToAminoAcid("Sumoylation1", AminoAcids.K);
		//Trimethylation	K
		addModificationToAminoAcid("Trimethylation", AminoAcids.K);
		//ubiquitinylation	K
		addModificationToAminoAcid("ubiquitinylation", AminoAcids.K);

	}
	
	public void addModificationToAminoAcid(String modificationName, byte aminoAcid) {
		ArrayList<Modification> modList =  getModificationList(aminoAcid);
		modList.add(getModification(modificationName));
	}
	
	@SuppressWarnings("unchecked")
	public ArrayList<Modification> getModificationList(byte aminoAcid) {
		return (ArrayList<Modification>) modificationsArray[aminoAcid];
	}
	
	public static Modification getModification(String modificationName) {
		for (Modification modification: Definitions.modifications) {
			if (modification.getPSI_MSname().equalsIgnoreCase(modificationName)) {
				return modification;
			}
		}
		U.p("No modification with this name found: " + modificationName);
		return null;
	}
	
	

}
