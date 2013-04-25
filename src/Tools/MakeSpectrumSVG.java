package Tools;

import java.io.File;
import java.util.ArrayList;

import Peppy.Match_IMP;
import Peppy.MatchesSpectrum;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.Spectrum;
import Peppy.SpectrumLoader;
import Peppy.U;
import Reports.MatchSVG;

public class MakeSpectrumSVG {

		public static void main(String [] args) {
			Properties.fragmentTolerance = 200;
			File spectrumFile = new File("/Users/risk2/PeppyData/ENCODE/GM12878/spectra uncompressed/wcl/15305/101111_GM_WCL_SDS_2_II.5801.5801.2.dta");
			Spectrum spectrum = SpectrumLoader.loadSpectra(spectrumFile).get(0);
			U.p("coverage: " + spectrum.getCoverage());
			U.p("spectrum mass: " + spectrum.getMass());
			Peptide peptide = new Peptide("APAGSAAGEGLLPHR");
			Match_IMP matchForReport = new Match_IMP();
			matchForReport.setPeptide(peptide);
			MatchesSpectrum matchesSpectrum = new MatchesSpectrum(spectrum);
			matchForReport.setSpectrumMatches(matchesSpectrum);
			ArrayList<Peppy.Match> matchesForReport = new ArrayList<Peppy.Match>();
			matchesForReport.add(matchForReport);
			
			File spectrumReportFolder = new File("");
			spectrumReportFolder.mkdir();

			
			/* save SVG */
			MatchSVG makeSVG = new MatchSVG(matchForReport, new File(spectrum.getMD5() + ".svg"));
			makeSVG.saveSVG();
			
			U.p("done");
		}
}
