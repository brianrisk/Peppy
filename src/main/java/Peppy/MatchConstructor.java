package Peppy;


/**
 * Copyright 2013, Brian Risk
 *
 * @author Brian Risk
 */
public class MatchConstructor {

    Class<? extends Match> match;

    @SuppressWarnings("unchecked")
    public MatchConstructor(String className) {
        try {
            match = (Class<? extends Match>) Class.forName(className);
        } catch (ClassNotFoundException e) {
            U.p("A match classname was not recognized");
            e.printStackTrace();
            System.exit(1);
        }
    }

    public Match createMatch(MatchesSpectrum matchesSpectrum, Peptide peptide) {
        return createMatch(matchesSpectrum, peptide, true);
    }

    public Match createMatch(MatchesSpectrum matchesSpectrum, Peptide peptide, boolean calculateScore) {
        Match out = null;
        try {
            out = match.newInstance();
            out.setMatchesSpectrum(matchesSpectrum);
            out.setPeptide(peptide);
            if (calculateScore) out.calculateScore();
        } catch (InstantiationException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IllegalAccessException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return out;
    }

}
