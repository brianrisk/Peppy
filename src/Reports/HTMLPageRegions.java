package Reports;

import java.io.File;
import java.util.ArrayList;

import Peppy.Region;
import Utilities.U;

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
		printTH("sequence");
		printTH("start");
		printTH("stop");
		printTH("-log10(P)");
		for(int i = 0; i < regions.size(); i++) {
			Region region = regions.get(i);
			
			/*make sub-report for this particular region */
			File regionsFile = new File(destinationFile.getParent(), "regions");
			regionsFile.mkdirs();
			File regionFile = new File(regionsFile, i + ".html");
			HTMLPageRegion regionReport = new HTMLPageRegion(regionFile, region);
			regionReport.makePage();
			
			
			printTR();
			printTD("" + i);
			printTD(region.getSequence().getSequenceFile().getName());
			printTD("<a href=\"regions/" + regionFile.getName() + "\">" + region.getStartLocation() + "</a>");
			
			/* UCSC link */
			String link = UCSC.getLink(region.getStartLocation(), region.getStopLocation(), region.getSequence());
			printTD(region.getStopLocation() + " (<a href=\"" +link + "\">UCSC</a>)");
			
			printTD("" + region.getPValue());
			
			
		}
		print("</table>");
		printFooter();

	}

}
