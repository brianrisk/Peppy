package Navigator;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import Peppy.Spectrum;
import Peppy.SpectrumLoader;
import Peppy.U;

/**
 * to be further developed!
 * @author Brian Risk
 *
 */
public class SpectrumReport {
	
	
	/**
	 * creates a text file that is a histogram of precursor distribution
	 */
	public static void spectrumReport() {
		
		ArrayList<Spectrum> spectra = SpectrumLoader.loadSpectra();
		U.p("number of spectra: "  + spectra.size());
		int binMax = 200;
		int binMaxMinusOne = binMax - 1;
		int [] specrumPeakDensityHistogram = new int[binMax];
		int index;
		for (Spectrum spectrum: spectra) {
//			index = spectrum.getPeakCount()/5;
			index = (int) spectrum.getMass() / 50;
			if (index > binMaxMinusOne) index = binMaxMinusOne;
			specrumPeakDensityHistogram[index]++;
		}
		/* create the output file */
		try {
			PrintWriter pw = new PrintWriter(new FileWriter("specrumPeakDensityHistogram.txt"));
			for (index = 0; index < binMax; index++) {
				pw.println(((index+ 1) * 50) + "\t" + specrumPeakDensityHistogram[index]);
			}
			pw.flush();
			pw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

}
