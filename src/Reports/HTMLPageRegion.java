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
		String regionName = region.getSequence().getSequenceFile().getName() + ": " + region.getStartLocation() + " to " + region.getStopLocation();
		printHeader(regionName, "<script src=\"http://peppyresearch.com/js/processing.js\"></script>");
		printP("Sequence: " + region.getSequence().getSequenceFile().getName());
		printP("region: " + region.getStartLocation() + " to " + region.getStopLocation());
		printP("region E value: " + region.getEValue());
		
		/* get access to the matches */
		ArrayList<Match> matches = region.getMatches();
		
		/* print out the processing section */
		int frameHeight = 15;
		int processngHeight = 6 * frameHeight;
		int processingWidth = 800;
		print("<script type=\"application/processing\">");
		print("void setup() {");
		print("noLoop();");
		print("size(" + processingWidth + ", " + processngHeight + ");");
		
		print("}");
		print("void draw() {");
		print("rect(0,0," + (processingWidth - 1) + " ," +  (processngHeight - 1) + ");");
		print("noStroke();");
		print("fill(0, 64);"); //a 25% black
		int x, y, width;
		y = 0;
		for(Match match: matches) {
			x = match.getPeptide().getStartIndex() - region.getStartLocation();
			x = scaleInt(x, region.getMaxLength(), processingWidth);
			y =  (match.getPeptide().getStartIndex() % 3) * frameHeight;
			if (match.getPeptide().isForward()) y += frameHeight * 3;
			width = match.getPeptide().getStopIndex() - match.getPeptide().getStartIndex();
			width = scaleInt(width, region.getMaxLength(), processingWidth);
			print("rect(" + x + ", " + y + ", " + width + ", " + frameHeight + ");");
		}
		print("}");
		print("</script><canvas width=\"" +  region.getMaxLength() + "px\" height=\"" + processngHeight + "px\"></canvas>");
		
		
		/* print out the table */
		print(ReportStrings.getTableHeader());
		
		NumberFormat nfDecimal = NumberFormat.getInstance();
		nfDecimal.setMaximumFractionDigits(4);
		NumberFormat nfPercent = NumberFormat.getPercentInstance();
		nfPercent.setMaximumFractionDigits(2);
		
		
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
	
	/**
	 * math func to let us scale coordinates when we are drawing a region smaller (or larger) than it is.
	 * @param x
	 * @param oldWidth
	 * @param newWidth
	 * @return
	 */
	private int scaleInt(int x, int oldWidth, int newWidth) {
		double dx = x;
		dx *= newWidth;
		dx /= oldWidth;
		return (int) Math.round(dx);
	}

}
