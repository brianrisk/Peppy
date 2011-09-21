package Reports;

import java.io.File;
import java.text.NumberFormat;
import java.util.ArrayList;

import Peppy.Match;
import Peppy.Properties;
import Peppy.Region;


public class HTMLPageRegion extends HTMLPage {
	
	Region region;
	
	public HTMLPageRegion(File destinationFile, Region region) {
		super(destinationFile);
		this.region = region;
	}

	@Override
	public void makePage() {
		printHeader();
		printP("Sequence: " + region.getSequence().getSequenceFile().getName());
		printP("region: " + region.getStartLocation() + " to " + region.getStopLocation());
		printP("region E value: " + region.getEValue());
		
		print(ReportStrings.getTableHeader());
		
		NumberFormat nfDecimal = NumberFormat.getInstance();
		nfDecimal.setMaximumFractionDigits(4);
		NumberFormat nfPercent = NumberFormat.getPercentInstance();
		nfPercent.setMaximumFractionDigits(2);
		
		ArrayList<Match> matches = region.getMatches();
		for(Match match: matches) {
			StringBuffer sb = new StringBuffer();
			sb.append("<tr>");
			
			sb.append("<td>");
			sb.append("<a href=\"spectra/");
			sb.append(match.getSpectrum().getId());
			sb.append(Properties.reportWebSuffix);
			sb.append("\">");
			sb.append(match.getSpectrum().getId());
			sb.append("</a>");
			sb.append("</td>");
			
			sb.append("<td>");
			sb.append(match.getPeptide().getAcidSequenceString());
			sb.append("</td>");
			
//			sb.append("<td>");
//			sb.append(match.getSpectrum().getFile().getAbsolutePath());
//			sb.append("</td>");
			
			
			sb.append("<td>");
			if (Properties.useSpliceVariants) {
				sb.append("null");
			} else {
				sb.append(match.getPeptide().getProtein().getName());
			}
			sb.append("</td>");
			
			sb.append("<td>");
			sb.append(match.getPeptide().getStartIndex());
			sb.append(" (" + match.getPeptide().getStartIndex() % 3 + ")");
			sb.append("</td>");
			
			sb.append("<td>");
			sb.append(match.getPeptide().getStopIndex());
			sb.append(" (" + match.getPeptide().getStopIndex() % 3 + ")");
			sb.append("</td>");
			
			sb.append("<td>");
			sb.append(match.getPeptide().isForward() ? "+" : "-");
			sb.append("</td>");
			
			sb.append("<td>");
			sb.append(match.getPeptide().isSpliced());
			sb.append("</td>");
			
			sb.append("<td>");
			sb.append(match.getSpectrum().getCharge());
			sb.append("</td>");
			
			sb.append("<td>");
			sb.append(match.getNumberOfModifications());
			sb.append("</td>");
			
			sb.append("<td>");
			sb.append(nfDecimal.format(match.getScore()));
			sb.append("</td>");
			
			sb.append("<td>");
			sb.append(match.getIonMatchTally());
			sb.append("</td>");
			
			sb.append("<td>");
			double matchPercent = (double) match.getIonMatchTally() / match.getPeptide().getAcidSequence().length;		
			sb.append(nfPercent.format(matchPercent));
			sb.append("</td>");
			
//			sb.append("<td>");
//			sb.append(nfDecimal.format(match.getScoreRatio()));
//			sb.append("</td>");
			
			sb.append("<td>");
			sb.append(match.getEValue());
			sb.append("</td>");
			
			//print out our table row
			print(sb.toString());
			
		}
		
		
		printFooter();

	}

}
