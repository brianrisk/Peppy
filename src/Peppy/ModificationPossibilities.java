package Peppy;

import java.util.ArrayList;

import Utilities.U;

public class ModificationPossibilities {
	
	ArrayList[] modificationsArray;
	
	public static void main(String args[]) {
		new ModificationPossibilities();
	}
	
	public ModificationPossibilities() {
		
		//initializes our modifications array
		modificationsArray = new ArrayList[AminoAcids.getNumberOfAminoAcids()];
		for (int i = 0; i < AminoAcids.getNumberOfAminoAcids(); i++) {
			ArrayList<Modification> mods = new ArrayList<Modification>();
			//add the fact that any amino acid may have NO modification
			mods.add(new Modification());
			modificationsArray[i] = mods;
		}
		
		//TODO remove later
		//specific to Yanbao

		//Acetylation	K, R
		addModificationToAminoAcid("Acetylation", AminoAcids.K);
		addModificationToAminoAcid("Acetylation", AminoAcids.R);
		//ADP ribosylation	E
		addModificationToAminoAcid("ADP Ribose addition", AminoAcids.E);
		//Biotinylation	K
		addModificationToAminoAcid("Biotinylation", AminoAcids.K);
		//Butyrylation	K
//		addModificationToAminoAcid("Butyrylation", AminoAcids.K);
		//Dimethylation	K, R
		addModificationToAminoAcid("Di-methylation", AminoAcids.K);
		addModificationToAminoAcid("Di-methylation", AminoAcids.R);
		//Methylation	K, R
		addModificationToAminoAcid("Methylation", AminoAcids.K);
		addModificationToAminoAcid("Methylation", AminoAcids.R);
		//Oxidation	M
		addModificationToAminoAcid("Oxidation or Hydroxylation", AminoAcids.M);
		//Palmitoylation	K
		addModificationToAminoAcid("Palmitoylation", AminoAcids.K);
		//Phosphorylation	S, T, Y
		addModificationToAminoAcid("Phosphorylation", AminoAcids.S);
		addModificationToAminoAcid("Phosphorylation", AminoAcids.T);
		addModificationToAminoAcid("Phosphorylation", AminoAcids.Y);
		//Propionylation	K
		addModificationToAminoAcid("Propionate labeling reagent light form (N-term & K)", AminoAcids.K);
		//Sumoylation1	K
		addModificationToAminoAcid("Sumo mutant Smt3-WT tail following trypsin digestion", AminoAcids.K);
		//Trimethylation	K
		addModificationToAminoAcid("Tri-methylation", AminoAcids.K);
		//ubiquitinylation	K
		addModificationToAminoAcid("ubiquitinylation residue", AminoAcids.K);

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
			if (modification.getDescription().equalsIgnoreCase(modificationName)) {
				return modification;
			}
		}
		U.p("No modification with this name found: " + modificationName);
		return null;
	}
	
	

}
