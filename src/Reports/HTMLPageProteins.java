package Reports;

import java.io.File;
import java.util.ArrayList;

import Peppy.Properties;
import Peppy.Protein;

public class HTMLPageProteins extends HTMLPage {
	
	private ArrayList<Protein> proteins;
	
	public HTMLPageProteins(ArrayList<Protein> proteins, File destinationFile) {
		super(destinationFile);
		this.proteins = proteins;
	}

	@Override
	public void makePage() {
		printHeader();
		
		printH1("Protein Report");
		printP("database: " + Properties.sequenceDirectoryOrFile);
		printP("spectra data set: " + Properties.spectraDirectoryOrFile);
		
		print("<table class=\"sortable\" id=\"box-table-a\">");
		print("<thead>");
		print("<tr>");
		printTH("name");
		printTH("coverage");
		printTH("area");
		printTH("mod coverage");
		printTH("mod area");
		print("</thead>");
		print("<tbody>");
		//make a directory to hold all of the proteins
		File proteinsDir = new File(destinationFile.getParentFile(), "proteins");
		for (int i = 0; i < proteins.size(); i++) {
			//make a page for our protein
			Protein protein = proteins.get(i);
			File proteinDir = new File(proteinsDir, protein.getName());
			proteinDir.mkdirs();
			File proteinFile = new File(proteinDir, "index.html");
			HTMLPageProtein php = new HTMLPageProtein(protein, proteinFile);
			php.makePage();
			
			//print out the table row
			//making sure that there is at least one match in the protein
			if (protein.getMatchesAll().size() > 0) {
				print("<tr>");
				printTD("<a href=\"proteins/" + protein.getName() + "/index.html\">" + protein.getName() + "</a>");
				printTD("" + protein.getMatchCoverage());
				printTD("" + protein.getMatchArea());
				printTD("" + protein.getModCoverage());
				printTD("" + protein.getModArea());
			}
		}
		print("</tbody>");
		print("</table>");
		
		
		printFooter();
	}

}
