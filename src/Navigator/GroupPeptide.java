package Navigator;

import java.util.ArrayList;

import Peppy.Location;
import Peppy.Match;
import Peppy.Peptide;
import Peppy.Spectrum;

public class GroupPeptide {
	
	Peptide peptide;
	
	/* all matches to this peptide */
	ArrayList<Match> matches;
	
	/* all spectra matching this peptide */
	ArrayList<Spectrum> spectra;
	
	/* all locations of this peptide */
	ArrayList<Location> locations;
	
	

}
