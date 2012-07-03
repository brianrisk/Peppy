package Reports;

import java.io.File;
import java.util.ArrayList;

import Peppy.Region;

/**
 * Copyright 2012, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class HTMLPageRegions extends HTMLPage {
	
	ArrayList<Region> regions;
	
	public HTMLPageRegions(File destinationFile, ArrayList<Region> regions) {
		super(destinationFile);
		this.regions = regions;
	}

	@Override
	public void makePage() {
		printHeader();
		print("<table class=\"sortable\" id=\"box-table-a\" width=\"95%\">");
		printTR();
		printTH("#");
		printTH("UCSC");
		printTH("sequence");
		printTH("start");
		printTH("hits");
		printTH("uniques");
		for(int i = 0; i < regions.size(); i++) {
			Region region = regions.get(i);
			
			/*make sub-report for this particular region */
			File regionsFile = new File(destinationFile.getParent(), "regions");
			regionsFile.mkdirs();
			File regionFile = new File(regionsFile, i + ".html");
			HTMLPageRegion regionReport = new HTMLPageRegion(regionFile, region);
			regionReport.makePage();
			
			
			printTR();
			
			/* region number */
			printTD("<a href=\"regions/" + regionFile.getName() + "\">" + i + "</a>");
			
			/* UCSC link */
			String link = UCSC.getLink(region.getStartLocation(), region.getStopLocation(), region.getSequence());
			printTD("(<a href=\"" +link + "\">UCSC</a>)");
			
			printTD(region.getSequence().getSequenceFile().getName());
			printTD("" + region.getStartLocation());

			
			printTD("" + region.getMatches().size());
			printTD("" + region.getUniqueCount());
			
			
		}
		print("</table>");
		printFooter();

	}

}
