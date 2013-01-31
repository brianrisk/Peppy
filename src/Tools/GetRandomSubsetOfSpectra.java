package Tools;

import java.io.File;
import java.util.ArrayList;
import java.util.Random;

import Peppy.Spectrum;
import Peppy.SpectrumLoader;
import Peppy.U;


/**
 * Let's say we have many, many spectra, but we want only a small subset to test.
 * 
 * This class lets you point at a directory of spectra and it will make a folder
 * that contains a random subset of all of the spectra in that directory.
 * 
 * @author Brian Risk
 *
 */
public class GetRandomSubsetOfSpectra {
	
	public static void main(String args[]) {
		/* setting up parameters */
//		File spectraFolder = new File("/Users/risk2/PeppyData/ENCODE/GM12878/spectra uncompressed/wcl/");
//		File destinationFolder = new File("UNC-ENCODE-WCL");
		
//		File spectraFolder = new File("/Users/risk2/PeppyData/UNC/spectra/CPTAC/UNC QExactive compRef/WHIM2");
//		File destinationFolder = new File("UNC-CPTAC-QEXATIVE");
		
		File spectraFolder = new File("/Users/risk2/PeppyData/WashU/spectra/43/WHIM2");
		File destinationFolder = new File("WASHU-CPTAC-TRIPLETOF");
		
		
		
		
		int numberOfSpectraForSubset = 2000;
		
		/* loading the spectra */
		U.p("loading spectra");
		ArrayList<Spectrum> spectra = SpectrumLoader.loadSpectraFromFolder(spectraFolder);
		
		/* get out if we don't have enough for a proper subset */
		if (spectra.size() <= numberOfSpectraForSubset) {
			U.p("There are only " + spectra.size() + " spectra to begin with!");
			return;
		}
		
		/* get the random subset */
		U.p("getting random subset");
		ArrayList<Spectrum> subset = new ArrayList<Spectrum>(numberOfSpectraForSubset);
		Random random = new Random();
		int randomIndex;
		for (int i = 0; i < numberOfSpectraForSubset; i++) {
			randomIndex = random.nextInt(spectra.size());
			subset.add(spectra.get(randomIndex));
			spectra.remove(randomIndex);
		}
		
		/* copy those spectra to the destination folder */
		U.p("copying spectra");
		destinationFolder.mkdirs();
		File destinationFile;
		for (Spectrum spectrum: subset) {
			destinationFile = new File(destinationFolder, spectrum.getFile().getName());
			U.copyfile(spectrum.getFile(), destinationFile);
		}
		
		U.p("done");
		
	}

}
