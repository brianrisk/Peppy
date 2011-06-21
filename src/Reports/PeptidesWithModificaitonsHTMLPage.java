package Reports;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;

import Peppy.Match;
import Peppy.Match_IMP_VariMod;
import Peppy.Properties;
import Peppy.Protein;

public class PeptidesWithModificaitonsHTMLPage extends HTMLPage {
	
	private ArrayList<Match_IMP_VariMod> matchesPTM;
	
	public PeptidesWithModificaitonsHTMLPage(ArrayList<Protein> proteins, File destinationFile) {
		super(destinationFile);
		matchesPTM = new ArrayList<Match_IMP_VariMod>();
		for (Protein protein: proteins) {
			matchesPTM.addAll(protein.getMatchesWithModifications());
		}
		Match.setSortParameter(Match_IMP_VariMod.SORT_BY_IMP_VALUE);
		Collections.sort(matchesPTM);
	}

	@Override
	public void makePage() {
		printHeader();
		
		printH1("Protein Report");
		printP("database: " + Properties.sequenceDirectoryOrFile);
		printP("spectra data set: " + Properties.spectraDirectoryOrFile);
		
		print("<table class=\"sortable\" id=\"box-table-a\" width=\"95%\">");
		print("<thead>");
		print("<tr>");
		printTH("#");
		printTH("peptide");
		printTH("protein");
		printTH("IMP");
		print("</thead>");
		print("<tbody>");
		
		//make a directory for the modifications
//		File modDir = new File(destinationFile.getParentFile(), "modifications");
//		modDir.mkdirs();
		
		Match_IMP_VariMod match;
		for (int i = 0; i < matchesPTM.size(); i++) {
			match = matchesPTM.get(i);
			
			//make a page for our match
//			File matchFile = new File(modDir, i + ".html");
//			MatchHTMLPage matchPage = new MatchHTMLPage(match, matchFile);
//			matchPage.makePage();
			
			//print out the table row
			print("<tr>");
			printTD("" + i);
			printTD("<a href=\"proteins/" + match.getPeptide().getProtein().getName() + "/" + match.getPeptide().getAcidSequenceString() + ".html\">" + match.getPeptide().getAcidSequenceString() + "</a>");
			printTD("" + match.getPeptide().getProtein().getName());
			printTD("" + match.getImpValue());
		}
		print("</tbody>");
		print("</table>");
		
		
		printFooter();
	}

}
