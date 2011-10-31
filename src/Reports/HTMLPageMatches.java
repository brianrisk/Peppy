package Reports;

import java.io.File;
import java.text.NumberFormat;
import java.util.ArrayList;

import Peppy.Match;
import Peppy.Match_IMP_VariMod;
import Peppy.Properties;
import Utilities.U;

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
		
		printMatchTable("spectra/");
		
		
		printFooter();
	}
	
	protected void printMatchTable(String spectraPath) {
		/* set up how we will round */
		NumberFormat nfDecimal = NumberFormat.getInstance();
		nfDecimal.setMaximumFractionDigits(2);
		
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
		printTH("e value");
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
		

		for (int i = 0; i < maxDisplay; i++) {
			/* get our match */
			Match match = matches.get(i);
			
			printTR();
			
			/* spectrum ID link */
			String spectrumLink = "<a href=\""  + spectraPath + match.getSpectrum().getId() + Properties.reportWebSuffix + "\">" + match.getSpectrum().getId() + "</a>";
			printTD(spectrumLink);
			
			/* the acid string */
			printTD(match.getPeptide().getAcidSequenceString());
			
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
			
			/* e value */
			printTD(nfDecimal.format(-Math.log10(match.getEValue())));
			
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
				printTD("" + match.hasMod());
			}
		}
		
		print("</table>");
	}

}
