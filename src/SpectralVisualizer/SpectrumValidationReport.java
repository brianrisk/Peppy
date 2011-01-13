package SpectralVisualizer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import Peppy.Peptide;
import Peppy.Spectrum;
import Utilities.U;

public class SpectrumValidationReport {
	
	public static void main(String args[]) {
		//load spectra
//		ArrayList<Spectrum> spectra = Spectrum.loadSpectraFromFolder("spectra encode membrane/GO_mem_FASP_dta201006");
		ArrayList<Spectrum> spectra = Spectrum.loadSpectraFromFolder("spectra encode membrane/SDS");
		
		//load list of spectra and peptides
		File inFile = new File("infilename.txt");
		ArrayList<String> spectraNames = new ArrayList<String>();
		ArrayList<String> peptideStrings = new ArrayList<String>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(inFile));
			String line = br.readLine();
			int lineNumber = 1;
			while (line != null) {
				String [] chunks = line.split("\t");
				if (chunks == null) {
					U.p("chunks is null");
					U.p("line number: " + lineNumber);
					break;
				}
				if (chunks.length != 2) {
					U.p("chunks length isn't 2");
					U.p("line number: " + lineNumber);
					break;
				}
				spectraNames.add(chunks[0]);
				peptideStrings.add(chunks[1]);
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//find the spectra from the list
		ArrayList<Spectrum> spectraSelect = new ArrayList<Spectrum>();
		ArrayList<Peptide> peptides = new ArrayList<Peptide>();
		for (int i = 0; i < spectraNames.size(); i++) {
			String name = spectraNames.get(i);
			for (Spectrum spectrum: spectra) {
				if (spectrum.getFile().getName().equals(name)) {
					spectraSelect.add(spectrum);
					peptides.add(new Peptide(peptideStrings.get(i)));
					break;
				}
			}
		}
		
		//make the report
		SpectralVisualizer.generateFullReport(spectraSelect, peptides);
	}

}
