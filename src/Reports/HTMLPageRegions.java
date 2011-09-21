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
		printP("<img src=\"histogram.jpg\">");
		print("<table class=\"sortable\" id=\"box-table-a\" width=\"95%\">");
		printTR();
		printTH("sequence");
		printTH("start");
		printTH("stop");
		printTH("E value");
		for(int i = 0; i < regions.size(); i++) {
			Region region = regions.get(i);
			
			/*make sub-report for this particular region */
			File regionsFile = new File(destinationFile.getParent(), "regions");
			regionsFile.mkdirs();
			File regionFile = new File(regionsFile, i + ".html");
			HTMLPageRegion regionReport = new HTMLPageRegion(regionFile, region);
			regionReport.makePage();
			
			
			printTR();
			printTD(region.getSequence().getSequenceFile().getName());
			printTD("<a href=\"regions/" + regionFile.getName() + "\">" + region.getStartLocation() + "</a>");
			
			/* UCSC link */
			String link = "http://genome.ucsc.edu/cgi-bin/hgTracks?hgHubConnect.destUrl=..%2Fcgi-bin%2FhgTracks&clade=mammal&org=Human&db=hg19&position=";
			link += U.getFileNameWithoutSuffix(region.getSequence().getSequenceFile());
			link += "%3A";
			link += region.getStartLocation();
			link += "-";
			link += region.getStopLocation();
			link += "&hgt.suggest=&hgt.suggestTrack=knownGene&&hgt.newJQuery=1&pix=922";
			printTD(region.getStopLocation() + " (<a href=\"" +link + "\">UCSC</a>)");
			
			printTD("" + region.getEValue());
			
			
		}
		print("</table>");
		printFooter();

	}

}
