package Reports;

import java.util.ArrayList;

import Peppy.Match;
import Peppy.Matches;
import Peppy.Properties;
import Peppy.Sequence;
import Peppy.U;


/** an object to help construct UCSC links
 * 
 * Copyright 2012, Brian Risk
 * Released under the Netscape Public License
 * 
 * @author Brian Risk
 *
 */
public class UCSC {
	
	/* okay, we'll maybe use this later when it is more general */
	private String clade; //eg. mammal
	private String organism; //eg. human
	private String database; //eg. hg19
	
	private static int speciesTracker = 0;
	public final static int HUMAN = speciesTracker++;
	public final static int MOUSE = speciesTracker++;
	
	public static String getLink(Match match) {
		return getLink(match.getPeptide().getStartIndex(), match.getPeptide().getStopIndex(), match.getPeptide().getParentSequence());
	}
	
	public static String getLink(int start, int stop, Sequence sequence) {
		if (Properties.isSequenceFileDNA) {
			String link = "http://genome.ucsc.edu/cgi-bin/hgTracks?hgHubConnect.destUrl=..%2Fcgi-bin%2FhgTracks&" + Properties.UCSCdatabase + "&position=";
			link += U.getFileNameWithoutSuffix(sequence.getSequenceFile());
			link += "%3A";
			link += start - 100;
			link += "-";
			link += stop + 100;
	//		link += "&hgt.suggest=&hgt.suggestTrack=knownGene&&hgt.newJQuery=1&pix=922";
			return link;
		} else {
			return "";
		}
	}
	
	public static String getLink(int start, int stop, String sequenceName) {
		String link = "http://genome.ucsc.edu/cgi-bin/hgTracks?hgHubConnect.destUrl=..%2Fcgi-bin%2FhgTracks&" + Properties.UCSCdatabase + "&position=";
		link += sequenceName;
		link += "%3A";
		link += start - 100;
		link += "-";
		link += stop + 100;
//		link += "&hgt.suggest=&hgt.suggestTrack=knownGene&&hgt.newJQuery=1&pix=922";
		return link;
	}
	
	public static String getButton(int start, int stop, Sequence sequence, ArrayList<Match> matches) {
		ArrayList<Match> bestMatches = Matches.getBestMatches(matches);
		String chromosome = U.getFileNameWithoutSuffix(sequence.getSequenceFile());
		StringBuffer out = new StringBuffer();
		out.append("<FORM ACTION=\"http://genome.ucsc.edu/cgi-bin/hgCustom\" METHOD=\"POST\"  ENCTYPE=\"multipart/form-data\" NAME=\"mainForm\" >");
		out.append("<INPUT TYPE=HIDDEN NAME='clade' VALUE='Mammal'>");
		out.append("<INPUT TYPE=HIDDEN NAME='org' VALUE='Human'>");
		out.append("<INPUT TYPE=HIDDEN NAME='db' VALUE='hg19'>");
		out.append("<INPUT TYPE=HIDDEN NAME='hgct_customText' VALUE=\"");
		out.append('\r');
		out.append("browser position " + chromosome + ":" + start + "-" + stop);
		out.append('\r');
		out.append("track name=PeppyRegion description='Peppy Region' visibility=2");
		out.append('\r');
		for (Match match: bestMatches) {
			out.append(chromosome + " " + match.getPeptide().getStartIndex() + " " + match.getPeptide().getStopIndex() + " " + match.getPeptide().getAcidSequenceString());
			out.append('\r');
		}
		out.append("\">");
		out.append('\r');
		out.append("<input type=\"submit\" name=\"Submit\" value=\"UCSC\" >");
		out.append("</form>");
		return out.toString();
	}

}
