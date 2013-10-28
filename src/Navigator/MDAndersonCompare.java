package Navigator;

import java.io.File;
import java.util.ArrayList;
import java.util.Hashtable;

import Peppy.Peppy;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.Protein;
import Peppy.Sequence_Protein;
import Peppy.Spectrum;
import Peppy.SpectrumLoader;
import Peppy.U;

/**
 * Delete this later.
 * 
 * My recent MD Anderson searches are not returning as many
 * as I did when I visited.  This is an attempt to find the
 * differences.
 * @author Brian Risk
 *
 */
public class MDAndersonCompare {
	
	public static void main(String args[]) {
		Peppy.init(args);
		findMissingMatch();
		U.p("done");
	}
	
	public static void saveSpecificSpectrum() {
		Properties.spectraDirectoryOrFile = new File("/Users/risk2/PeppyData/Hanash/spectra/HoustonApril/IPAS_0506_CID_Surface_SG31to32.mgf");
		ArrayList<Spectrum> spectra = SpectrumLoader.loadSpectra();
		for(Spectrum spectrum: spectra) {
			if (spectrum.getMD5().equals("c04bac9bf554fb63eb37a0f49918cd87")) {
				spectrum.saveDTA(new File("/Users/risk2/Documents/workspace/JavaGFS/"));
				break;
			}
		}
	}
	
	public static void checkThatPeptideIsCreated() {
		String inQuestion = "LENLDSDVVQLR";
		
		Sequence_Protein db = new Sequence_Protein(new File("/Users/risk2/PeppyData/public/sequences/protein/UniProt_Human_2012_03.fasta")); 
		ArrayList<Peptide> peptides = db.extractAllPeptides(false);
		for (Peptide peptide: peptides) {
			if (peptide.getAcidSequenceString().equals(inQuestion)) {
				U.p("found " + inQuestion + "!");
				break;
			}
		}
	}
	
	public static void simpleBLAST() {
		String proteinSequence = "MAACTARRALAVGSRWWSRSLTGARWPRPLCAAAGAGAFSPASTTTTRRHLSSRNRPEGKVLETVGVFEVPKQNGKYETGQLFLHSIFGYRGVVLFPWQARLYDRDVASAAPEKAENPAGHGSKEVKGKTHTYYQVLIDARDCPHISQRSQTEAVTFLANHDDSRALYAIPGLDYVSHEDILPYTSTDQVPIQHELFERFLLYDQTKAPPFVARETLRAWQEKNHPWLELSDVHRETTENIRVTVIPFYMGMREAQNSHVYWWRYCIRLENLDSDVVQLRERHWRIFSLSGTLETVRGRGVVGREPVLSKEQPAFQYSSHVSLQASSGHMWGTFRFERPDGSHFDVRIPPFSLESNKDEKTPPSGLHW";
		Protein protein = new Protein("Q9Y2S7", proteinSequence, false);
		Peptide peptide = new Peptide("LENLDSDVVQLR");
	}
	
	
	public static void findMissingMatch() {
		File resultsPast = new File("/Users/risk2/Documents/MD Anderson/Visit results/group1/1 group1 - UniProt_Human_2012_03.fasta/report.txt");
		File resultsPresent = new File("/Users/risk2/Documents/workspace/JavaGFS/reports/houstonApril/1 IPAS_0506_CID_Surface_SG31to32.mgf - UniProt_Human_2012_03.fasta/report.txt");
		
		ArrayList<Match> matchesPast = Match.loadMatches(resultsPast);
		ArrayList<Match> matchesPresent = Match.loadMatches(resultsPresent);
		
		//reduce the past machtes to only the fraction we are examining
		ArrayList<Match> matchesPastReduced = new ArrayList<Match>();
		for(Match match: matchesPast) {
			File spectrumFile = match.getFile("FilePath");
			if (spectrumFile.getName().endsWith("Surface_SG31to32.mgf")) {
				matchesPastReduced.add(match);
			}
		}
		matchesPast = matchesPastReduced;
		
		//create hash
		Hashtable<String, Match> matchesPresentHash = new Hashtable<String, Match> (); 
		for(Match match: matchesPresent) {
			String peptide = match.getString("peptideSequence");
			if (matchesPresentHash.get(peptide) == null) {
				matchesPresentHash.put(peptide,  match);
			}
		}
		
		//go though what we had and report the first one from the present that doesn't exist
		for (Match matchPast: matchesPast) {
			String peptidePast = matchPast.getString("peptideSequence");
			if (matchesPresentHash.get(peptidePast) == null) {
				
				//exceptions
				if (peptidePast.equals("LENLDSDVVQLR")) continue;
				if (peptidePast.equals("LAGEELAGEEAPQEK")) continue;
				if (peptidePast.equals("VPGISSIEQGMTGLK")) continue;
				if (peptidePast.equals("LENLDSDVVQLR")) continue;
				
				
				// print report
				U.p("Match past:");
				U.p("peptide: " + peptidePast);
				U.p("score: " + matchPast.getDouble("score"));
				break;
			}
		}
	}

}
