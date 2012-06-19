package Tools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import Peppy.Peak;
import Peppy.Spectrum;
import Peppy.SpectrumLoader;
import Peppy.U;

/**
 * How many format conversion tools will I write in my lifetime?
 * Here I create a tool to take the Spectrum object and save it as
 * a MGF file for Mascot searches
 * 
 * Why am I even making this an object and not a static method?
 * I don't know.  It's early.  Just... just go over there.
 * @author Brian Risk
 *
 */
public class SpectrumToMGF {
	
	private File outputFile;
	private Spectrum spectrum;
	
	/**
	 * This main is for the specific jobs I am doing
	 * @param args
	 */
	public static void main(String args[]) {
		//read in our array of spectra
		String testDirectoryName = "/Users/risk2/PeppyData/tests/";
		
		//all you need to change is this one testName variable
		String testName = "human";
		
		
		ArrayList<Spectrum> spectra = SpectrumLoader.loadSpectraFromFolder(testDirectoryName + testName + "/spectra");
		
		//save properties
		File saveFolder = new File ("MGF/" + testName);
		saveFolder.mkdirs();
		
		//for each, create SpectrumToMGF object, output
		for (Spectrum spectrum: spectra) {
			File saveLocation = new File(saveFolder, spectrum.getFile().getName() + ".mgf");
			SpectrumToMGF mgf = new SpectrumToMGF(spectrum, saveLocation);
			mgf.generateMGF();
		}
		
		U.p("done");
	}
	
	public SpectrumToMGF(Spectrum spectrum, File outputFile) {
		this.spectrum = spectrum;
		this.outputFile = outputFile;
	}
	
	public void generateMGF() {
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));
			ArrayList<Peak> peaks = spectrum.getPeaks();
			pw.println("BEGIN IONS");
			pw.println("TITLE=" + spectrum.getFile().getName());
			pw.println("PEPMASS=" + spectrum.getMass());
			for (Peak peak: peaks) {
				pw.println(peak.getMass() + " " + peak.getIntensity());
			}
			pw.println("END IONS");
			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
