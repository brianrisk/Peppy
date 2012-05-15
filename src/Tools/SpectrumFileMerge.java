package Tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import Peppy.U;

/**
 * Built because some search tools use one file containing all spectra rather than many spectrum files
 * @author Brian Risk
 *
 */
public class SpectrumFileMerge {
	
	//takes the spectra from a folder and merges them into a file
	public static void main(String args[]) {
		
		File destination;
		File folder;
		try {
			U.p("starting");
			
			//Here I was merging the original spectra files
//			//E. coli
//			destination = new File("ecoli.pkl");
//			folder = new File("/Users/risk2/PeppyData/tests/ecoli/spectra");
//			merge(folder, destination, ".txt");
//			
//			//Kapp
//			destination = new File("kapp.dta");
//			folder = new File("/Users/risk2/PeppyData/tests/human/spectra");
//			merge(folder, destination, ".dta");
//			
//			//Aurum
//			destination = new File("aurum.pkl");
//			folder = new File("/Users/risk2/PeppyData/tests/aurum/spectra");
//			merge(folder, destination, ".pkl");
			
			//Here I am merging the MGF files I made
			//E. coli
			destination = new File("ecoli.mgf");
			folder = new File("MGF/ecoli");
			merge(folder, destination, ".mgf");
			
			//Kapp
			destination = new File("kapp.mgf");
			folder = new File("MGF/human");
			merge(folder, destination, "mgf");
			
			//Aurum
			destination = new File("aurum.mgf");
			folder = new File("MGF/aurum");
			merge(folder, destination, ".mgf");
			
			U.p("done");
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void merge(File spectrumFolder, File destination, String extension) throws IOException {
		//set up our destination file
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(destination)));
		//read through all the files and concatenate them
		File [] spectra = spectrumFolder.listFiles();
		String line;
		for (int i = 0; i < spectra.length; i++) {
			if (!spectra[i].getName().toLowerCase().endsWith(extension)) continue;
			BufferedReader br = new BufferedReader(new FileReader(spectra[i]));
			line = br.readLine();
			while (line != null) {
				line = line.trim();
				if (line != "") {
					pw.println(line);
				}
				line = br.readLine();
			}
			pw.println();
		}
		pw.flush();
		pw.close();
	}

}
