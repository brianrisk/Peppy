package Reports;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;

import Peppy.Match;
import Peppy.Properties;
import Peppy.Sequence;
import Peppy.Spectrum;
import SpectralVisualizer.SpectralVisualizer;
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
	 * @param sequences
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
					if (Peppy.Properties.useEValueCutOff) {
						if (match.getEValue() < Peppy.Properties.eValueCutOff) {
							bestMatches.add(match);
						}
					} else {
						bestMatches.add(match);
					}
				} else {
					matchRank++;
				}
//				match.setRank(matchRank);
				if (matchRank == 2) {
					scoreRatio = previousMatch.getScore() / match.getScore();
					previousMatch.setTandemFitScoreRatio(scoreRatio);
				}
				previousMatch = match;
			}
			
			//sort our best matches by score ratio
//			Match.setSortParameter(Match.SORT_BY_SCORE_RATIO);
			Match.setSortParameter(Match.SORT_BY_E_VALUE);
			Collections.sort(bestMatches);
			
			//set up our main index file
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(indexFile)));
			U.appendFile(pw, Properties.reportWebHeaderFile);
			pw.println("<h1>Best match for each spectrum</h1>");
			pw.println("<p>Number of spectra with quality matches: " + bestMatches.size() + "</p>");
			U.appendFile(pw, Properties.reportWebTableHeader);
			
			
			for (int i = 0; i < bestMatches.size(); i++) {
				Match match = bestMatches.get(i);
				StringBuffer sb = new StringBuffer();
				sb.append("<tr>");
				
				//create animation
//				U.p(i);
//				new TandemFitAnimation(match.getSpectrum().getFile().getAbsolutePath(), match.getPeptide().getAcidSequence());
				
				sb.append("<td>");
				sb.append("<a href=\"spectra/");
				sb.append(match.getSpectrum().getId());
				sb.append(Properties.reportWebSuffix);
				sb.append("\">");
				sb.append(match.getSpectrum().getId());
				sb.append("</a>");
				sb.append("</td>");
				
				sb.append("<td>");
				sb.append(match.getPeptide().getAcidSequence());
				sb.append("</td>");
				
//				sb.append("<td>");
//				sb.append(match.getSpectrum().getFile().getAbsolutePath());
//				sb.append("</td>");
				
				sb.append("<td><nobr>");
				sb.append("<a href=\"sequences/");
				sb.append(match.getSequence().getId());
				sb.append(Properties.reportWebSuffix);
				sb.append("\">");
				sb.append(match.getSequence().getSequenceFile().getName());
				sb.append("</a> ");
				sb.append("</nobr></td>");
				
				sb.append("<td>");
				sb.append("<a href=\"neighborhoods/");
				sb.append(i);
				sb.append(Properties.reportWebSuffix);
				sb.append("\">");
				sb.append(match.getPeptide().getStartIndex());
				sb.append(" (" + match.getPeptide().getStartIndex() % 3 + ")");
				sb.append("</a>");
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
				sb.append((double) match.getIonMatchTally() / match.getPeptide().getAcidSequence().length());
				sb.append("</td>");
				
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
			
			U.appendFile(pw, Properties.reportWebFooterFile);
			pw.flush();
			pw.close();
			
			if (Peppy.Properties.isSequenceFileDNA &&  Properties.generateSequenceReport) {
				for (int i = 0; i < sequences.size(); i++) {
					generateSequenceReport(sequences.get(i));
				}
			}
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
			if (thisMatch.getSequence() != match.getSequence()) continue;
			if (Math.abs(thisMatch.getPeptide().getStartIndex() - match.getPeptide().getStartIndex()) > Properties.locusNeighborhood) continue;
			theseMatches.add(thisMatch);
		}

		try {
			neighborhoodDirectory.mkdirs();
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(indexFile)));
			U.appendFile(pw, Properties.reportWebHeaderSubFile);
			
			//print the best match for each spectrum
			pw.println("<h1>Matches in the neighborhood of " + match.getPeptide().getStartIndex() + "</h1>");
			U.appendFile(pw, Properties.reportWebTableHeader);
			Collections.sort(theseMatches);
			for (int i = 0; i < theseMatches.size(); i++) {
				pw.println(getTableRow(theseMatches.get(i), i));
			}
			
			U.appendFile(pw, Properties.reportWebFooterFile);
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
	
	public void generateSequenceReport(Sequence sequence) {
		/*
		 * The main index page will have a list with each spectrum and their top match.
		 * If there are more than one sequence then we will list each sequence and 
		 * how many matches for each spectra there were.
		 */
		
		
		File sequenceReportDirectory = new File(reportDir, "sequences");
		sequenceReportDirectory.mkdirs();
		File indexFile = new File(sequenceReportDirectory, sequence.getId() + Properties.reportWebSuffix);
		ArrayList<Match> theseMatches = getMatchesWithSequence(sequence, matches);
		try {
			
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(indexFile)));
			U.appendFile(pw, Properties.reportWebHeaderSubFile);
			
			//region report
//			insertSequenceRegionReport(theseMatches, sequenceReportDirectory, pw, sequence);
			GeneReport geneReport = new GeneReport(sequence, theseMatches, 3000);
			//draw the e value histogram
			pw.println("<h2>E value histogram for genes</h2>");
			pw.println("<p>");
			File histogramFile = new File(sequenceReportDirectory, sequence.getId() + "-hist.jpg");
			HistogramVisualizer.drawHistogram(geneReport.getEValueHistogram(), 300, 300, histogramFile);
			pw.println("<img src=\"" + histogramFile.getName() + "\">");
			pw.println("</p>");
			geneReport.insertSequenceRegionReport(theseMatches, sequenceReportDirectory, pw);
			
			//print the best match for each spectrum
//			//collect the best matches
//			ArrayList<Match> bestMatches = new ArrayList<Match>();
//			for (int i = 0; i < spectra.size(); i++) {
//				ArrayList<Match> specific = getMatchesWithSpectrum(spectra.get(i), theseMatches);
//				Collections.sort(specific);
//				if (specific.size() > 0)
//					bestMatches.add(specific.get(0));
//			}
//			Collections.sort(bestMatches);
//			pw.println("<h1>Best match for each spectrum in " + sequence.getSequenceFile().getName() + "</h1>");
//			appendFile(pw, Properties.reportWebTableHeader);
//			
//			for (int i = 0; i < bestMatches.size(); i++) {
//				pw.println(getTableRow(bestMatches.get(i), i));
//			}
			
			U.appendFile(pw, Properties.reportWebFooterFile);
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
	
	/**
	 * Give this method a print writer and some other necessary ingredients and it will put the report in there
	 * 
	 * 
	 * @param theseMatches
	 * @param sequenceReportDirectory
	 * @param pw
	 * @param sequence
	 * @throws IOException 
	 */
//	private void insertSequenceRegionReport(ArrayList<Match> theseMatches, File sequenceReportDirectory, PrintWriter pw, Sequence sequence) throws IOException {
//		//initialize our variables
//		int numberOfTopRegions = 10;
//		int regionSize = 3000;
//		int histogramHeight = 200;
//		int sequenceLength = sequence.getSequenceLength();
//		int numberOfRegions = (int) Math.ceil(sequenceLength / regionSize) + 1;
//		ArrayList<SequenceRegion> regions = new ArrayList<SequenceRegion>(numberOfRegions);
//		for (int i = 0 ; i < numberOfRegions; i++) {
//			regions.add(new SequenceRegion(regionSize * i, regionSize));
//		}
//		int matchIndex;
//		SequenceRegion targetRegion;
//		
//		//fill our regions with matches
//		for (Match match: theseMatches) {
//			matchIndex = match.getPeptide().getIndex();
//			targetRegion = regions.get(matchIndex / regionSize);
//			targetRegion.addHit(match);
//		}
//		Collections.sort(regions);
//		int maxBar = 0;
//		SequenceRegion region;
//		for (int i = 0; i < numberOfTopRegions; i++) {
//			region = regions.get(i);
//			if (region.getMaxBar() > maxBar) maxBar = region.getMaxBar();
//		}
//		
//		//create visualizations and report of the top regions
//		File regionDirectory = new File(sequenceReportDirectory, "regions" + sequence.getId());
//		regionDirectory.mkdirs();
//		pw.println("<table border=1>");
//		for (int i = 0; i < numberOfTopRegions; i++) {
//			region = regions.get(i);
//			
//			pw.println();
//			pw.print("<tr>");
//			
//			pw.print("<td>");
//			pw.print("<a href=\"" + regionDirectory.getName() + "/" + i + ".html\">");
//			pw.print(region.getStartIndex());
//			pw.print("</a>");
//			pw.print("</td>");
//			
//			pw.print("<td>");
//			pw.print(region.getScore());
//			pw.print("</td>");
//			
//			pw.print("<td>");
//			pw.print(region.getMaxBar());
//			pw.print("</td>");
//			
//			pw.print("<td>");
//			pw.print("<img src=\"" + regionDirectory.getName() + "/" + i + ".jpg\">");
//			pw.print("</td>");
//			
//			pw.print("</tr>");
//			
//			//create the images and html for this region
//			HistogramVisualizer.drawShadedHistogram(region.getHistogram(), histogramHeight, maxBar, new File(regionDirectory, i + ".jpg"));
//		}
//		pw.println("</table>");
//	}
	
	public void generateSpectrumReport(Spectrum spectrum) {
		File sequenceDirectory = new File(reportDir, "spectra");
		File indexFile = new File(sequenceDirectory, spectrum.getId() + Properties.reportWebSuffix);
		ArrayList<Match> theseMatches = getMatchesWithSpectrum(spectrum, matches);
		try {
			sequenceDirectory.mkdirs();
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(indexFile)));
			U.appendFile(pw, Properties.reportWebHeaderSubFile);
			
			pw.println("<h1>" + spectrum.getFile().getName() + "</h1>");
			
			//draw the e value histogram
			pw.println("<h2>E value histogram for spectrum " + spectrum.getId() + "</h2>");
			pw.println("<p>");
			File histogramFile = new File(sequenceDirectory, spectrum.getId() + "-hist.jpg");
			HistogramVisualizer.drawHistogram(spectrum.getEValueCalculator().getHistogram(), 300, 300, histogramFile);
			pw.println("<img src=\"" + histogramFile.getName() + "\">");
			pw.println("</p>");
			
			//sorting
			Match.setSortParameter(Match.SORT_BY_SCORE);
			Collections.sort(theseMatches);
			
			//draw spectrum visualizations for top 2 matches
			pw.println("<h2>Ion matches for top two spectra</h2>");
			pw.println("<p>");
			int atLeastTwo = 2;
			if (theseMatches.size() < atLeastTwo) atLeastTwo = theseMatches.size();
			for (int i = 0; i < atLeastTwo; i++) {
				Match match = theseMatches.get(i);
				File spectrumVisualization = new File(sequenceDirectory, spectrum.getId() + "-spect-" + i + ".jpg");
				SpectralVisualizer.markMatchingIons(spectrum, match.getPeptide());
				SpectralVisualizer.drawSpectrum(spectrum, 500, 200, spectrumVisualization, true);
				pw.println("<a href=\"" + spectrumVisualization.getName() + "\"><img src=\"" + spectrumVisualization.getName() + "\" border=0></a><br>");
			}
			pw.println("</p>");
			
			//print the best match for each spectrum
			pw.println("<h2>Matches for spectrum " + spectrum.getId() + "</h2>");
			U.appendFile(pw, Properties.reportWebTableHeader);
			for (int i = 0; i < theseMatches.size(); i++) {
				pw.println(getTableRow(theseMatches.get(i), i));
			}
			
			U.appendFile(pw, Properties.reportWebFooterFile);
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
	
	private ArrayList<Match> getMatchesWithSpectrum(Spectrum spectrum, ArrayList<Match> theseMatches) {
		ArrayList<Match> out = new ArrayList<Match>();
		for (int i = 0; i < theseMatches.size(); i++) {
			Match match = theseMatches.get(i);
			if (match.getSpectrum().getId() == spectrum.getId()) {
				out.add(match);
			}
		}
		return out;
	}
	
	private ArrayList<Match> getMatchesWithSequence(Sequence sequence, ArrayList<Match> theseMatches) {
		ArrayList<Match> out = new ArrayList<Match>();
		for (int i = 0; i < theseMatches.size(); i++) {
			Match match = theseMatches.get(i);
			if (match.getSequence().getId() == sequence.getId()) {
				out.add(match);
			}
		}
		return out;
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
		sb.append(match.getPeptide().getAcidSequence());
		sb.append("</td>");
		
		sb.append("<td><nobr>");
		if (Properties.generateSequenceReport) {
			sb.append("<a href=\"../sequences/");
			sb.append(match.getSequence().getId());
			sb.append(Properties.reportWebSuffix);
			sb.append("\">");
		}
		sb.append(match.getSequence().getSequenceFile().getName ());
		if (Properties.generateSequenceReport) {
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
		sb.append((double) match.getIonMatchTally() / match.getPeptide().getAcidSequence().length());
		sb.append("</td>");
		
		sb.append("<td>");
		sb.append(match.getEValue());
		sb.append("</td>");
		
		return sb.toString();
	}

}
