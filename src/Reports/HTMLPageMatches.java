package Reports;

import java.io.File;
import java.text.NumberFormat;
import java.util.ArrayList;

import Peppy.Match;
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
		
		printMatchTable();
		
		
		printFooter();
	}
	
	protected void printMatchTable() {
		/* set up how we will round */
		NumberFormat nfDecimal = NumberFormat.getInstance();
		nfDecimal.setMaximumFractionDigits(2);
		
		/* print our table headers */
		print("<table class=\"sortable\" id=\"box-table-a\" width=\"95%\">");
		printTR();
		printTH("ID");
		printTH("UCSC");
		printTH("acid");
		if (!Properties.useSpliceVariants) {
			printTH("sequence");
		}
		printTH("start");
		printTH("stop");
		printTH("strand");
		if (Properties.useSpliceVariants) {
			printTH("spliced");
		}
		printTH("charge");
		printTH("score");
		printTH("e value");
		printTH("Hydrophobic percentage");
		printTH("Hydrophilic percentage");
		if (Properties.useIsotopeLabeling) {
			printTH("isotope confirmed");
		}
		
		for (int i = 0; i < maxDisplay; i++) {
			/* get our match */
			Match match = matches.get(i);
			
			printTR();
			
			/* spectrum ID link */
			String spectrumLink = "<a href=\"spectra/" + match.getSpectrum().getId() + Properties.reportWebSuffix + "\">" + match.getSpectrum().getId() + "</a>";
			printTD(spectrumLink);
			
			/* UCSC link */
			String link = UCSC.getLink(match);
			printTD("(<a href=\"" +link + "\">UCSC</a>)");
			
			/* the acid string */
			printTD(match.getPeptide().getAcidSequenceString());
			
			/* the sequence / protein name */
			if (!Properties.useSpliceVariants) {
				printTD(match.getPeptide().getProtein().getName());
			}
			
			/* start */
			printTD("" + match.getPeptide().getStartIndex());
			
			/* stop */
			printTD("" + match.getPeptide().getStopIndex());
			
			/* strand */
			printTD(match.getPeptide().isForward() ? "+" : "-");
			
			/* is spliced? */
			if (Properties.useSpliceVariants) {
				printTD("" + match.getPeptide().isSpliced());
			}
			
			/* spectrum charge */
			printTD("" + match.getSpectrum().getCharge());
			
			/* score */
			printTD(nfDecimal.format(match.getScore()));
			
			/* e value */
			printTD(nfDecimal.format(-Math.log10(match.getEValue())));
			
			/* hydrophobic */
			printTD(nfDecimal.format(match.getPeptide().getHydrophobicProportion()));
			
			/* hydrophilic */
			printTD(nfDecimal.format(match.getPeptide().getHydrophilicProportion()));
			
			if (Properties.useIsotopeLabeling) {
				printTD("" + match.isHasIsotopeConfirmation());
			}
		}
		
		print("</table>");
	}

}
