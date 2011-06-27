package Reports;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;

import Peppy.Match;
import Peppy.Properties;
import Peppy.Sequence;
import Peppy.Spectrum;
import Utilities.U;


/**
 * Okay, so you've got all of your results.  Now what?
 * I'll tell you now what.  You want to see them presented in
 * a nice, easy and intuitive manner.  That's what this
 * class does.
 * @author Brian Risk
 *
 */
public class HTMLReporter {
	
	ArrayList<Match> matches;
	ArrayList<Spectrum> spectra;
	ArrayList<Sequence> sequences;
	File reportDir;
	
	
	/**
	 * @param matches
	 * @param spectra
	 * @param sequence_DNAs
	 */
	public HTMLReporter(ArrayList<Match> matches,
			ArrayList<Spectrum> spectra, ArrayList<Sequence> sequences, File reportDir) {
		this.matches = matches;
		this.spectra = spectra;
		this.sequences = sequences;
		this.reportDir = reportDir;
	}


	public void generateFullReport() {
		//create our report directory
		reportDir.mkdirs();
		File indexFile = new File(reportDir, "index" + Properties.reportWebSuffix);
		try {
			
			
			//sorting our matches by spectrum then score
			Match.setSortParameter(Match.SORT_BY_SPECTRUM_ID_THEN_SCORE);
			Collections.sort(matches);
			
			/* 
			 * finding best matches
			 * finding ratio of best match in a spectrum to the runner-up
			 * setting    rank
			 * set HMM Score
			 */
			ArrayList<Match> bestMatches = new ArrayList<Match>();
			

			int spectrumID = -1; //an ID will never be -1
			int matchRank = 0;
			double scoreRatio;
			Match previousMatch = null;
			for (Match match: matches) {
				if (spectrumID != match.getSpectrum().getId()) {
					spectrumID = match.getSpectrum().getId();
					matchRank = 1;
					if (match.getEValue() < Peppy.Properties.maxEValue) {
						bestMatches.add(match);
					}
				} else {
					matchRank++;
				}
//				match.setRank(matchRank);
				if (matchRank == 2) {
					scoreRatio = previousMatch.getScore() / match.getScore();
					previousMatch.setScoreRatio(scoreRatio);
				}
				previousMatch = match;
			}
			
			//sort our best matches by score ratio
//			Match.setSortParameter(Match.SORT_BY_SCORE_RATIO);
			Match.setSortParameter(Match.SORT_BY_IMP_VALUE);
			Collections.sort(bestMatches);
			
			//set up our main index file
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(indexFile)));
			//HTML header
			pw.println(ReportStrings.getHeader());
			
			//Title
			pw.println("<h1>Best match for each spectrum</h1>");
			pw.println("<p>Number of spectra with quality matches: " + bestMatches.size() + "</p>");
			
			//table header
			pw.println(ReportStrings.getTableHeader());
			
			NumberFormat nfDecimal = NumberFormat.getInstance();
			nfDecimal.setMaximumFractionDigits(4);
			NumberFormat nfPercent = NumberFormat.getPercentInstance();
			nfPercent.setMaximumFractionDigits(2);
			
			//limit how many we will display
			int maxDisplay = 10000;
			if (maxDisplay > bestMatches.size()) maxDisplay = bestMatches.size();
			
			for (int i = 0; i < maxDisplay; i++) {
				Match match = bestMatches.get(i);
				StringBuffer sb = new StringBuffer();
				sb.append("<tr>");
				
				sb.append("<td>");
				sb.append("<a href=\"spectra/");
				sb.append(match.getSpectrum().getId());
				sb.append(Properties.reportWebSuffix);
				sb.append("\">");
				sb.append(match.getSpectrum().getId());
				sb.append("</a>");
				sb.append("</td>");
				
				sb.append("<td>");
				sb.append(match.getPeptide().getAcidSequenceString());
				sb.append("</td>");
				
//				sb.append("<td>");
//				sb.append(match.getSpectrum().getFile().getAbsolutePath());
//				sb.append("</td>");
				
				
				sb.append("<td>");
				sb.append(match.getPeptide().getProtein().getName());
				sb.append("</td>");
				
				sb.append("<td>");
//				sb.append("<a href=\"neighborhoods/");
//				sb.append(i);
//				sb.append(Properties.reportWebSuffix);
//				sb.append("\">");
				sb.append(match.getPeptide().getStartIndex());
				sb.append(" (" + match.getPeptide().getStartIndex() % 3 + ")");
//				sb.append("</a>");
				sb.append("</td>");
				
				sb.append("<td>");
				sb.append(match.getPeptide().getStopIndex());
				sb.append(" (" + match.getPeptide().getStopIndex() % 3 + ")");
				sb.append("</td>");
				
				sb.append("<td>");
				sb.append(match.getPeptide().isForward() ? "+" : "-");
				sb.append("</td>");
				
				sb.append("<td>");
				sb.append(match.getPeptide().isSpliced());
				sb.append("</td>");
				
				sb.append("<td>");
				sb.append(match.getSpectrum().getCharge());
				sb.append("</td>");
				
				sb.append("<td>");
				sb.append(match.getNumberOfModifications());
				sb.append("</td>");
				
				sb.append("<td>");
				sb.append(nfDecimal.format(match.getScore()));
				sb.append("</td>");
				
				sb.append("<td>");
				sb.append(match.getIonMatchTally());
				sb.append("</td>");
				
				sb.append("<td>");
				double matchPercent = (double) match.getIonMatchTally() / match.getPeptide().getAcidSequence().length;		
				sb.append(nfPercent.format(matchPercent));
				sb.append("</td>");
				
//				sb.append("<td>");
//				sb.append(nfDecimal.format(match.getScoreRatio()));
//				sb.append("</td>");
				
				sb.append("<td>");
				sb.append(match.getEValue());
				sb.append("</td>");
				
				//print out our table row
				pw.println(sb);
				
				//generate neighborhood reports
				if (Peppy.Properties.isSequenceFileDNA && Properties.generateNeighborhoodReport) {
					generateNeighborhoodReport(bestMatches.get(i), i);
				}
			}
			
			pw.println(ReportStrings.getFooter());
			pw.flush();
			pw.close();
			
			if (Properties.generateSpectrumReport) {
				for (Match match: bestMatches) {
					generateSpectrumReport(match.getSpectrum());
				}
			}
			
		} catch (FileNotFoundException e) {
			U.p("could not find file: " + indexFile.getName());
			e.printStackTrace();
		} catch (IOException e) {
			U.p("could not read file: " + indexFile.getName());
			e.printStackTrace();
		}
	}
	
	public void generateNeighborhoodReport(Match match, int index) {
		File neighborhoodDirectory = new File(reportDir, "neighborhoods");
		File indexFile = new File(neighborhoodDirectory, index + Properties.reportWebSuffix);
		
		//find all matches that are from the same chromosome
		ArrayList<Match> theseMatches = new ArrayList<Match>();
		for (int i = 0; i < matches.size(); i++) {
			Match thisMatch = matches.get(i);
			if (thisMatch.getPeptide().getParentSequence() != match.getPeptide().getParentSequence()) continue;
			if (Math.abs(thisMatch.getPeptide().getStartIndex() - match.getPeptide().getStartIndex()) > Properties.locusNeighborhood) continue;
			theseMatches.add(thisMatch);
		}

		try {
			neighborhoodDirectory.mkdirs();
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(indexFile)));
			pw.println(ReportStrings.getHeader());
			
			//print the best match for each spectrum
			pw.println("<h1>Matches in the neighborhood of " + match.getPeptide().getStartIndex() + "</h1>");
			pw.println(ReportStrings.getTableHeader());
			Collections.sort(theseMatches);
			for (int i = 0; i < theseMatches.size(); i++) {
				pw.println(getTableRow(theseMatches.get(i), i));
			}
			
			pw.println(ReportStrings.getFooter());
			pw.flush();
			pw.close();
			
		} catch (FileNotFoundException e) {
			U.p("could not find file: " + indexFile.getName());
			e.printStackTrace();
		} catch (IOException e) {
			U.p("could not read file: " + indexFile.getName());
			e.printStackTrace();
		}
	}
	
	
	public void generateSpectrumReport(Spectrum spectrum) {
		//set up our files and folders
		File sequenceDirectory = new File(reportDir, "spectra");
		sequenceDirectory.mkdirs();
		File indexFile = new File(sequenceDirectory, spectrum.getId() + Properties.reportWebSuffix);
		
		//our report object
		HTMLPage page;
		
		//select which type of spectrum report to give
		if (Properties.scoringMethodName.equals("Peppy.Match_IMP_MultiMod")) {
			page = new SpectrumMultiMod(spectrum, matches, indexFile);
		} else {
			page = new SpectrumHTMLPage(spectrum, matches, indexFile);
		}
		page.makePage();
	}
	
	

	private String getTableRow(Match match, int index) {
		StringBuffer sb = new StringBuffer();
		sb.append("<tr>");
		
		sb.append("<td>");
		if (Properties.generateSpectrumReport) {
			sb.append("<a href=\"../spectra/");
			sb.append(match.getSpectrum().getId());
			sb.append(Properties.reportWebSuffix);
			sb.append("\">");
		}
		sb.append(match.getSpectrum().getId());
		if (Properties.generateSpectrumReport) {
			sb.append("</a>");
		}
		sb.append("</td>");
		
		sb.append("<td>");
		sb.append(match.getPeptide().getAcidSequenceString());
		sb.append("</td>");
		
		sb.append("<td><nobr>");
		if (Properties.generateSequenceReport) {
			sb.append("<a href=\"../sequences/");
			sb.append(match.getPeptide().getParentSequence().getId());
			sb.append(Properties.reportWebSuffix);
			sb.append("\">");
			sb.append(match.getPeptide().getParentSequence().getSequenceFile().getName ());
			sb.append("</a>");
		}
		if (match.getPeptide().isForward()) {sb.append(" forward ");}
		else {sb.append(" reverse ");}
		sb.append("</nobr></td>");
		
		
		sb.append("<td>");
		sb.append(match.getPeptide().getStartIndex());
		sb.append(" (" + match.getPeptide().getStartIndex() % 3 + ")");
		sb.append("</td>");
		
		sb.append("<td>");
		sb.append(match.getPeptide().getStopIndex());
		sb.append(" (" + match.getPeptide().getStopIndex() % 3 + ")");
		sb.append("</td>");
		
		sb.append("<td>");
		sb.append(match.getPeptide().isForward() ? "+" : "-");
		sb.append("</td>");
		
		sb.append("<td>");
		sb.append(match.getPeptide().isSpliced());
		sb.append("</td>");
		
		sb.append("<td>");
		sb.append(match.getScore());
		sb.append("</td>");
		
		sb.append("<td>");
		sb.append(match.getIonMatchTally());
		sb.append("</td>");
		
		sb.append("<td>");
		sb.append((double) match.getIonMatchTally() / match.getPeptide().getAcidSequence().length);
		sb.append("</td>");
		
		sb.append("<td>");
		sb.append(match.getEValue());
		sb.append("</td>");
		
		return sb.toString();
	}

}
