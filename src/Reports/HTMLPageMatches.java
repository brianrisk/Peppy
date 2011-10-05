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
		
		for (int i = 0; i < maxDisplay; i++) {
			/* get our match */
			Match match = matches.get(i);
			
			printTR();
			
			/* spectrum ID link */
			String spectrumLink = "<a href=\"spectra/" + match.getSpectrum().getId() + Properties.reportWebSuffix + "\">" + match.getSpectrum().getId() + "</a>";
			printTD(spectrumLink);
			
			/* UCSC link */
			//String link = "http://genome.ucsc.edu/cgi-bin/hgTracks?hgHubConnect.destUrl=..%2Fcgi-bin%2FhgTracks&clade=mammal&org=Mouse&db=mm9&position=";
			String link = "http://genome.ucsc.edu/cgi-bin/hgTracks?hgHubConnect.destUrl=..%2Fcgi-bin%2FhgTracks&clade=mammal&org=Human&db=mm9&position=";
			link += U.getFileNameWithoutSuffix(match.getPeptide().getParentSequence().getSequenceFile());
			link += "%3A";
			link += match.getPeptide().getStartIndex();
			link += "-";
			link += match.getPeptide().getStartIndex();
			link += "&hgt.suggest=&hgt.suggestTrack=knownGene&&hgt.newJQuery=1&pix=922";
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
		}
		
		print("</table>");
	}

}
