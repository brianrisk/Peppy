package Experimental;

import Navigator.MatchRow;
import Peppy.*;
import Reports.MatchSVG;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;

public class SpectrumDeNovo {
	
	Spectrum spectrum;
	ArrayList<Peak> peaks;
	ArrayList<Double> masses;
	protected String sequence;
	
	double deltaMin = 57.01; // 57.02146
	double deltaMax = 186.09; //186.07931
	double fragmentTolerance = .03;
	

	/* testing it out */
	public static void main(String args[]) {
		U.p("loading matches");
		ArrayList<MatchRow> matches = new ArrayList<MatchRow>();
//		matches.addAll(MatchRow.loadMatches(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM16-CompRef/1 WHIM16 - gencode proteome/report.txt")));
//		matches.addAll(MatchRow.loadMatches(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM16-CompRef/4 WHIM16 - mouse proteome/report.txt")));
//		matches.addAll(MatchRow.loadMatches(new File("/Users/risk2/PeppyData/CPTAC/reports/UNC-WHIM16-CompRef/6 WHIM16 - varimod/report.txt")));
		
		
//		matches.addAll(MatchRow.loadMatches(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/UNC-ENCODE-WCL/1 UNC-ENCODE-WCL - UniProt_Human_2012_03.fasta/report.txt")));
//		matches.addAll(MatchRow.loadMatches(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/UNC-CPTAC-QEXATIVE/1 UNC-CPTAC-QEXATIVE - UniProt_Human_2012_03.fasta/report.txt")));
		
		matches.addAll(MatchRow.loadMatches(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/WASHU-CPTAC-TTOF/1 WASHU-CPTAC-TTOF - UniProt_Human_2012_03.fasta/report.txt")));
//		matches.addAll(MatchRow.loadMatches(new File("/Users/risk2/Documents/workspace/JavaGFS/reports/KAPP/1 spectra - UniProt_Human_2012_03.fasta/report.txt")));
		Collections.sort(matches);
		
		MatchRow match = matches.get(0);
		ArrayList<Spectrum> spectra = SpectrumLoader.loadSpectraFromFile(match.getFile("FilePath"));
		SpectrumDeNovo spectrum = new SpectrumDeNovo(spectra.get(0));
	
		U.p(match.getString("peptideSequence"));
		U.p(spectrum.deNovo());
//		spectrum.deNovo();
		
		Match_Blank matchBlank = new Match_Blank(spectra.get(0), new Peptide(match.getString("peptideSequence")), match.getScore());
		MatchSVG spectrumVisualization = new MatchSVG(matchBlank, new File ("fancy psm.svg"));
		spectrumVisualization.saveSVG();
	}
	
	
	public SpectrumDeNovo(Spectrum spectrum) {
		this.spectrum = spectrum;	
		peaks = spectrum.getPeaks();
		
		masses = new ArrayList<Double>(peaks.size());
		for (Peak peak: peaks) {
			masses.add(new Double(peak.getMass()));
		}
		
		removeCounterpartIons();
	}
	
	public String deNovo() {
		ArrayList<StringBuffer> strings = new ArrayList<StringBuffer>();
		for (int indexA = 0; indexA < peaks.size() - 1; indexA ++) {
			if (peaks.get(indexA).used == false) {
				ArrayList<StringBuffer> suffixStrings = recursiveDeNovo(indexA, 0);
				strings.addAll(suffixStrings);
				
			}
		}
		
		/* find the longest string length */
		int longestStringLength = 0;
		for (StringBuffer sb: strings) {
			if (sb.length() > longestStringLength) longestStringLength = sb.length();
		}
		
		/* return the fist longest string */
		for (StringBuffer sb: strings) {
			if (sb.length() == longestStringLength) {
				U.p(sb);
				return sb.toString();
			}
		}
		return null;
	}
	
	
	/**
	 * Uses recursion to find de novo sequence
	 * @param peakIndex
	 * @return
	 */
	private ArrayList<StringBuffer> recursiveDeNovo(int peakIndex, double aggregateError) {
		double massA, massB, delta, error;
		ArrayList<StringBuffer> strings = new ArrayList<StringBuffer>();
		
		massA = peaks.get(peakIndex).getMass();
		for (int indexB = peakIndex + 1; indexB < peaks.size(); indexB ++) {
			massB = peaks.get(indexB).getMass();
			delta = massB - massA;
			if (delta > deltaMax) break;
			if (delta < deltaMin) continue;
			Character acid =  getAcid(delta);
			if (acid != '*') {
				error = (delta - AminoAcids.getWeightMono(acid)) + aggregateError;
				if (Math.abs(aggregateError) > fragmentTolerance) continue;
				ArrayList<StringBuffer> suffixStrings = recursiveDeNovo(indexB, error);
				if (suffixStrings.size() != 0) {
					for (StringBuffer suffix: suffixStrings) {
						StringBuffer add = new StringBuffer();
						add.append(acid);
						add.append(suffix);
						strings.add(add);
						peaks.get(indexB).used = true;
					}
				} else {
					StringBuffer add = new StringBuffer();
					add.append(acid);
					strings.add(add);
				}
			}
		}
	
		
		/* find the longest string length */
		int longestStringLength = 0;
		for (StringBuffer sb: strings) {
			if (sb.length() > longestStringLength) longestStringLength = sb.length();
		}
		
		/* count how many have the longest length */
		int longStringCount = 0;
		for (StringBuffer sb: strings) {
			if (sb.length() == longestStringLength) longStringCount++;
		}
		
		/* make a new array with this amount */
		ArrayList<StringBuffer> out = new ArrayList<StringBuffer>(longStringCount);
		
		/* populate the array with the strings of this length */
		for (StringBuffer sb: strings) {
			if (sb.length() == longestStringLength) out.add(sb);
//			U.p(sb);
		}
		
		/* return the reduced array */
		return out;
			
	}
	
	private void removeCounterpartIons() {
		int foundIndex;
		double counterpart, foundMass;
//		U.p(peaks.size());
		double massStop = spectrum.getMass() / 2;
		for (double mass: masses) {
			if (mass > massStop) break;
			counterpart = getCounterpartMass(mass);
			foundIndex = Collections.binarySearch(masses, counterpart);
			if (foundIndex >=0 ) {
				/* the odd chance that the counter part mass is *exactly* matched */
				peaks.remove(foundIndex);
			} else {
				foundIndex *= -1;
				foundIndex -= 1;
				
				if (foundIndex < 0) continue;
				if (foundIndex >= masses.size()) continue;
				
				U.p(foundIndex + ", " + peaks.size());
				foundMass = masses.get(foundIndex);
				if (Math.abs(foundMass - counterpart) < fragmentTolerance) {
					peaks.remove(foundIndex);
					continue;
				}
				
				foundIndex += 2;
				if (foundIndex < 0) continue;
				if (foundIndex >=  masses.size()) continue;
				foundMass = masses.get(foundIndex);
				if (Math.abs(foundMass - counterpart) < fragmentTolerance) {
					peaks.remove(foundIndex);
					continue;
				}
			}
		}
//		U.p(peaks.size());
		
	}
	
	private char getAcid(double mass) {
		int roundedMass = (int) Math.round(mass);
		if (roundedMass == 57) {if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.G)) < fragmentTolerance) return 'G';}
		if (roundedMass == 71) {if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.A)) < fragmentTolerance) return 'A';}
		if (roundedMass == 87) {if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.S)) < fragmentTolerance) return 'S';}
		if (roundedMass == 97) {if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.P)) < fragmentTolerance) return 'P';}
		if (roundedMass == 99) {if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.V)) < fragmentTolerance) return 'V';}
		if (roundedMass == 101) {if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.T)) < fragmentTolerance) return 'T';}
		if (roundedMass == 103) {if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.C)) < fragmentTolerance) return 'C';}
		if (roundedMass == 113) {if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.L)) < fragmentTolerance) return 'L';} /* or I */
		if (roundedMass == 114) {if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.N)) < fragmentTolerance) return 'N';}
		if (roundedMass == 115) {if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.D)) < fragmentTolerance) return 'D';}
		if (roundedMass == 128) {if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.K)) < fragmentTolerance) return 'K';}
		if (roundedMass == 128) {if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.Q)) < fragmentTolerance) return 'Q';}
		if (roundedMass == 129) {if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.E)) < fragmentTolerance) return 'E';}
		if (roundedMass == 131) {if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.M)) < fragmentTolerance) return 'M';}
		if (roundedMass == 137) {if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.H)) < fragmentTolerance) return 'H';}
		if (roundedMass == 147) {if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.F)) < fragmentTolerance) return 'F';}
		if (roundedMass == 156) {if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.R)) < fragmentTolerance) return 'R';}
		if (roundedMass == 163) {if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.Y)) < fragmentTolerance) return 'Y';}
		if (roundedMass == 186) {if (Math.abs(mass - AminoAcids.getWeightMono(AminoAcids.W)) < fragmentTolerance) return 'W';}
		return '*';
	}


	private double getCounterpartMass(double mass) {
		//NOTE:  would these two hydrogens be for +2 charge or what? 
		return spectrum.getMass() - mass + Definitions.HYDROGEN_MONO + Definitions.HYDROGEN_MONO;
	}
	

}
