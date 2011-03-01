package Reports;

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;

import Peppy.Match;
import Peppy.Protein;

public class ProteinHTMLPage extends HTMLPage{
	
	static Color PINK = new Color(255, 128, 128);
	private Protein protein;
	
	public ProteinHTMLPage(Protein protein, File destinationFile) {
		super(destinationFile);
		this.protein = protein;
	}
	
	@Override
	public void makePage() {
		printHeader();
		
		//printing basic info
		String acidString = protein.getAcidString();
		int [] matchPositions = protein.getMatchPositions();
		printH3(protein.getName());
		if (matchPositions == null) {
			printP("Nothing to report");
		} else {
			printP("Match coverage: " + protein.getMatchCoverage());
			printP("Protein score: " + protein.getScore());
			
			
			//printing the acid table which shows match regions
			print("<table width=\"80%\">");
			for (int i = 0; i < acidString.length(); i++) {
				if (i % 20 == 0) {
					print("<tr>");
					printTD("" + i);
				}
				if (matchPositions[i] == Protein.T_FPR01) printTD("" + acidString.charAt(i), Color.RED);
				if (matchPositions[i] == Protein.T_FPR05) printTD("" + acidString.charAt(i), PINK);
				if (matchPositions[i] == Protein.T_FPRXX) printTD("" + acidString.charAt(i), Color.LIGHT_GRAY);
				if (matchPositions[i] == Protein.T_MOD) printTD("" + acidString.charAt(i), Color.GREEN);
				if (matchPositions[i] == Protein.T_NOTHING) printTD("" + acidString.charAt(i));
			}
			print("</table>");
			
			//printing the table of matches
			ArrayList<Match> matches = protein.getMatchesAll();
			Match.setSortParameter(Match.SORT_BY_IMP_VALUE);
			Collections.sort(matches);
			print("<table class=\"sortable\" id=\"box-table-a\">");
			print("<thead>");
			print("<tr>");
			printTH("acid sequence");
			printTH("start");
			printTH("IMP");
			printTH("score");
			printTH("has mod");
			print("</thead>");
			print("<tbody>");
			Match match;
			for (int i = 0; i < matches.size(); i++) {
				match = matches.get(i);
				if (match.hasModification()) {
					print("<tr bgcolor=\"#00FF00\">");
				} else {
					print("<tr>");	
				}
				
				printTD("<a href=\"" + match.getPeptide().getAcidSequenceString() + ".html\">" + match.getPeptide().getAcidSequenceString() + "</a>");
				printTD("" + match.getPeptide().getStartIndex());
				printTD("" + match.getImpValue());
				printTD("" + match.getScore());
				printTD("" + match.hasModification());
				
				//make a page for the match
				File matchFile = new File(destinationFile.getParentFile(), match.getPeptide().getAcidSequenceString() + ".html");
				MatchHTMLPage matchPage = new MatchHTMLPage(match, matchFile);
				matchPage.makePage();
			}
			print("</tbody>");
			print("</table>");
		}
		printFooter();
	}
	

}
