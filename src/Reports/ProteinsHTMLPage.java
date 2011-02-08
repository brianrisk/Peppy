package Reports;

import java.io.File;
import java.util.ArrayList;

import Peppy.Properties;
import Peppy.Protein;

public class ProteinsHTMLPage extends HTMLPage {
	
	private ArrayList<Protein> proteins;
	
	public ProteinsHTMLPage(ArrayList<Protein> proteins, File destinationFile) {
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
		for (int i = 0; i < proteins.size(); i++) {
			//make a page for our protein
			Protein protein = proteins.get(i);
			File proteinDir = new File(destinationFile.getParentFile(), "protein" + i);
			File proteinFile = new File(proteinDir, "index.html");
			proteinDir.mkdirs();
			ProteinHTMLPage php = new ProteinHTMLPage(protein, proteinFile);
			php.makePage();
			
			//print out the table row
			print("<tr>");
			printTD("<a href=\"protein" + i + "/index.html\">" + protein.getName() + "</a>");
			printTD("" + protein.getMatchCoverage());
			printTD("" + protein.getMatchArea());
			printTD("" + protein.getModCoverage());
			printTD("" + protein.getModArea());
		}
		print("</tbody>");
		print("</table>");
		
		
		printFooter();
	}

}
