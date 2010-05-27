package SpectrumComparison;

import Peppy.Peptide;
import Peppy.Spectrum;

/**
 * A data structure to hold both a spectrum and the peptide sequence which
 * produced the spectrum.
 * @author Brian Risk
 *
 */
public class SpectrumPeptidePair {
	
	Spectrum spectrum;
	Peptide peptide;
	
	public SpectrumPeptidePair(Spectrum spectrum, Peptide peptide) {
		this.spectrum = spectrum;
		this.peptide = peptide;
	}

	public Spectrum getSpectrum() {
		return spectrum;
	}

	public void setSpectrum(Spectrum spectrum) {
		this.spectrum = spectrum;
	}

	public Peptide getPeptide() {
		return peptide;
	}

	public void setPeptide(Peptide sequence) {
		this.peptide = sequence;
	}
	
	

}
