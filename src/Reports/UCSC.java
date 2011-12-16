package Reports;

import Peppy.Match;
import Peppy.Properties;
import Peppy.Sequence;
import Utilities.U;


/** an object to help construct UCSC links
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
			String link = "http://genome.ucsc.edu/cgi-bin/hgTracks?hgHubConnect.destUrl=..%2Fcgi-bin%2FhgTracks&clade=mammal&org=Human&db=hg19&position=";
			if (Properties.isYale) {
				link = "http://genome.ucsc.edu/cgi-bin/hgTracks?hgHubConnect.destUrl=..%2Fcgi-bin%2FhgTracks&clade=mammal&org=Mouse&db=mm9&position=";
			}
			link += U.getFileNameWithoutSuffix(sequence.getSequenceFile());
			link += "%3A";
			link += start;
			link += "-";
			link += stop;
	//		link += "&hgt.suggest=&hgt.suggestTrack=knownGene&&hgt.newJQuery=1&pix=922";
			return link;
		} else {
			return "";
		}
	}

}
