package Peppy;


/**
 * Every Spectrum object will be associated with a collection of matches.
 * <p>
 * This class will be the judge and jury for which matches become associated
 * with a spectrum.  Things such as confidence thresholds, eliminating duplicates
 * and whatnot.
 * <p>
 * Thing become so very simple when we say, hey, we want only the best matches
 * for every spectrum.  Then we don't have to worry about ranks and rank counts
 * and removing duplicate matches or all of the terrible memory issues.
 * <p>
 * Further, we have the added benefit of not needing a minimum score as every
 * spectrum will have
 * <p>
 * Copyright 2013, Brian Risk
 *
 * @author Brian Risk
 */
public class MatchesSpectrum extends Matches {

    Spectrum spectrum;


    public MatchesSpectrum(Spectrum spectrum) {
        this.spectrum = spectrum;
    }

    public Spectrum getSpectrum() {
        return spectrum;
    }


}
