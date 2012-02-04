package Reports;

import java.io.File;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;

import Peppy.Match;
import Peppy.Matches;
import Peppy.Peptide;
import Peppy.Properties;

public class HTMLPageMatches extends HTMLPage {
	
	ArrayList<Match> matches;
	int maxDisplay;

	public HTMLPageMatches(ArrayList<Match> matches, int maxDisplay, File destinationFile) {
		super(destinationFile);
		this.matches = matches;
		this.maxDisplay = maxDisplay;
	}

	@Override
	public void makePage() {
		
		printHeader();
		
		printLink("regions html/index.html", "View regions report");
		
		/* print table of best matches */
		printH2("Best matches for each peptide");
		printBestMatchTable("spectra/");
		
		printFooter();
	}
	
	protected void printBestMatchTable(String spectraPath) {
		printMatchTable(spectraPath, Matches.getBestMatches(matches));
	}
	
	protected void printMatchTable(String spectraPath) {
		printMatchTable(spectraPath, matches);
	}
	
	protected void printMatchTable(String spectraPath, ArrayList<Match> matches) {
		/* print our table headers */
		print("<table class=\"sortable\" id=\"box-table-a\" width=\"95%\">");
		printTR();
		printTH("ID");
		printTH("acid");
		printTH("location");
		if (Properties.useSpliceVariants) {
			printTH("spliced");
		}
		printTH("score");
		if (Properties.isYale) {
			printTH("inORF");
			printTH("Hydrophobic percentage");
			printTH("Hydrophilic percentage");
		}
		if (Properties.useIsotopeLabeling) {
			printTH("isotope confirmed");
		}
		
		if (Properties.searchModifications) {
			printTH("has mod");
		}
		printTH("rank");
		printTH("count");
		
		/*print the rows */
		printMatchRows(matches, spectraPath);	
		
		print("</table>");
	}
	
	
	protected void printMatchRows(ArrayList<Match> matches, String spectraPath) {
		/* some basic bounds checks */
		if (matches == null) return;
		if (matches.size() == 0) return;
		int maxDisplay = this.maxDisplay;
		if (matches.size() < maxDisplay) maxDisplay = matches.size();
		
		/* set up how we will round */
		NumberFormat nfDecimal = NumberFormat.getInstance();
		nfDecimal.setMaximumFractionDigits(2);
		
		for (int i = 0; i < maxDisplay; i++) {
			/* get our match */
			Match match = matches.get(i);
			
			printTR();
			
			/* spectrum ID link */
			String spectrumLink = "<a href=\""  + spectraPath + match.getSpectrum().getId() + Properties.reportWebSuffix + "\">" + match.getSpectrum().getId() + "</a>";
			printTD(spectrumLink);
			
			/* the acid string */
			String peptideMarkup = "";
			if (match.rank > 1) peptideMarkup += "<font color=#444444>";
			peptideMarkup += match.getPeptide().getAcidSequenceString();
			if (match.rank > 1) peptideMarkup += "</font>";
			/* google link */
			String link = "http://www.google.com/search?hl=en&q=" + match.getPeptide().getAcidSequenceString();
			printTD("(<a href=\"" +link + "\">" + peptideMarkup + "</a>)");
			
			/* the sequence / protein name */
			if (!Properties.useSpliceVariants) {
				printTD(match.getPeptide().getProtein().getName() + " " + match.getPeptide().getStartIndex() + " " + ( match.getPeptide().isForward() ? "+" : "-"));
			}

			
			/* is spliced? */
			if (Properties.useSpliceVariants) {
				printTD("" + match.getPeptide().isSpliced());
			}
			
			/* score */
			printTD(nfDecimal.format(match.getScore()));
			
			if (Properties.isYale) {
				/* ORF */
				printTD("" + match.getPeptide().isInORF());
				
				/* hydrophobic */
				printTD(nfDecimal.format(match.getPeptide().getHydrophobicProportion()));
				
				/* hydrophilic */
				printTD(nfDecimal.format(match.getPeptide().getHydrophilicProportion()));
			}
			
			if (Properties.useIsotopeLabeling) {
				printTD("" + match.isHasIsotopeConfirmation());
			}
			
			/* PTM */
			if (Properties.searchModifications) {
				printTD("" + match.hasModification());
			}
			
			/* rank */
			printTD("" + match.rank);
			printTD("" + match.rankCount);
		}
	}

}
