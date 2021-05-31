package Reports;

import Peppy.Spectrum;
import Peppy.SpectrumLoader;
import Peppy.U;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

/**
 * to be further developed!
 * @author Brian Risk
 *
 */
public class SpectrumReport {
	
	public static void main(String [] argw) {
		U.p("creating spectrum statistics report");
		U.startStopwatch();
		Peppy.Peppy.init(argw);
		ArrayList<Spectrum> spectra = SpectrumLoader.loadSpectra();
		U.p("number of spectra: "  + spectra.size());
		spectrumReport(spectra, new File("spetrum statistics"));
		U.stopStopwatch();
	}
	
	
	/**
	 * creates a text file that is a histogram of precursor distribution
	 */
	public static void spectrumReport(ArrayList<Spectrum> spectra, File mainReportDir) {	
		File spectrumStatistics = new File(mainReportDir, "spectrum statistics");
		spectrumStatistics.mkdirs();
		
		int binMax = 200;
		int binMaxMinusOne = binMax - 1;
		int [] spectrumPeakDensityHistogram = new int[binMax];
		int [] spectrumMassDistribution = new int[binMax];
		int index;
		for (Spectrum spectrum: spectra) {
			index = spectrum.getPeakCount()/5;
			if (index > binMaxMinusOne) index = binMaxMinusOne;
			spectrumPeakDensityHistogram[index]++;
			
			index = (int) spectrum.getMass() / 50;
			if (index > binMaxMinusOne) index = binMaxMinusOne;
			spectrumMassDistribution[index]++;
		}
		/* create the output file */
		try {
			PrintWriter spectrumPeakDensityHistogramFile = new PrintWriter(new FileWriter(new File(spectrumStatistics, "specrumPeakDensityHistogram.txt")));
			for (index = 0; index < binMax; index++) {
				spectrumPeakDensityHistogramFile.println(((index+ 1) * 5) + "\t" + spectrumPeakDensityHistogram[index]);
			}
			spectrumPeakDensityHistogramFile.flush();
			spectrumPeakDensityHistogramFile.close();
			
			PrintWriter spectrumMassDistributionFile = new PrintWriter(new FileWriter(new File(spectrumStatistics, "spectrumMassDistribution.txt")));
			for (index = 0; index < binMax; index++) {
				spectrumMassDistributionFile.println(((index+ 1) * 50) + "\t" + spectrumMassDistribution[index]);
			}
			spectrumMassDistributionFile.flush();
			spectrumMassDistributionFile.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

}
