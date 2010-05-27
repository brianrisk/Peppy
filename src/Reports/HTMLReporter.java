package Reports;

import java.io.*;
import java.util.*;
import Peppy.*;
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
	
	ArrayList<SpectrumPeptideMatch> matches;
	ArrayList<Spectrum> spectra;
	ArrayList<Sequence> sequences;
	
	
	/**
	 * @param matches
	 * @param spectra
	 * @param sequences
	 */
	public HTMLReporter(ArrayList<SpectrumPeptideMatch> matches,
			ArrayList<Spectrum> spectra, ArrayList<Sequence> sequences) {
		this.matches = matches;
		this.spectra = spectra;
		this.sequences = sequences;
	}


	public void generateFullReport() {
		File indexFile = new File(Properties.reportDirectory, "index" + Properties.reportWebSuffix);
		try {
			//create our report directory
			Properties.reportDirectory.mkdirs();
			
			//set up our main index file
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(indexFile)));
			appendFile(pw, Properties.reportWebHeaderFile);
			pw.println("<h1>Best match for each spectrum</h1>");
			appendFile(pw, Properties.reportWebTableHeader);
			
			//sorting our matches by spectrum then score
			SpectrumPeptideMatch.setSortParameter(SpectrumPeptideMatch.SORT_BY_SPECTRUM_ID);
			Collections.sort(matches);
			
			/*
			 * finding best matches
			 * finding ratio of best match in a spectrum to the runner-up
			 * setting    rank
			 * set HMM Score
			 */
			ArrayList<SpectrumPeptideMatch> bestMatches = new ArrayList<SpectrumPeptideMatch>(spectra.size());
			
			SpectrumPeptideMatch match = matches.get(0);
			int matchRank = 1;
			int spectrumID = match.getSpectrum().getId();
			if (matches.size() > 1) {
				match.setMSMSFitScoreRatio(match.getScoreMSMSFit() / matches.get(1).getScoreMSMSFit());
			}
			match.setMSMSFitRank(matchRank);
			matchRank++;
			bestMatches.add(match);
			
			for (int i = 1; i < matches.size(); i++) {
				match = matches.get(i);
				if (match.getSpectrum().getId() != spectrumID) {
					bestMatches.add(match);
					spectrumID = matches.get(i).getSpectrum().getId();
					matchRank = 1;
					if (i + 1 < matches.size()) {
						match.setMSMSFitScoreRatio( match.getScoreMSMSFit() / matches.get(i + 1).getScoreMSMSFit());
					}
				}
				match.setMSMSFitRank(matchRank);
				matchRank++;
			}
			
			//sort our best matches by score ratio
			SpectrumPeptideMatch.setSortParameter(SpectrumPeptideMatch.SORT_BY_SCORE_RATIO);
			Collections.sort(bestMatches);
			
			for (int i = 0; i < bestMatches.size(); i++) {
				match = bestMatches.get(i);
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
				sb.append(match.getPeptide().getSequence());
				sb.append("</td>");
				
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
				sb.append(match.getPeptide().getIndex());
				sb.append("</a>");
				sb.append("</td>");
				
				sb.append("<td>");
				sb.append(match.getScoreHMM());
				sb.append("</td>");
				
				sb.append("<td>");
				sb.append(match.getMSMSFitScoreRatio());
				sb.append("</td>");
				
				sb.append("<td>");
				sb.append(match.getEValue());
				sb.append("</td>");
				
				//print out our table row
				pw.println(sb);
				
				//generate neighborhood reports
//				generateNeighborhoodReport(bestMatches.get(i), i);
			}
			
			appendFile(pw, Properties.reportWebFooterFile);
			pw.flush();
			pw.close();
			
//			for (int i = 0; i < sequences.size(); i++) {
//				generateSequenceReport(sequences.get(i));
//			}
//			
//			for (int i = 0; i < spectra.size(); i++) {
//				generateSpectrumReport(spectra.get(i));
//			}
			
		} catch (FileNotFoundException e) {
			U.p("could not find file: " + indexFile.getName());
			e.printStackTrace();
		} catch (IOException e) {
			U.p("could not read file: " + indexFile.getName());
			e.printStackTrace();
		}
	}
	
	public void generateNeighborhoodReport(SpectrumPeptideMatch match, int index) {
		File neighborhoodDirectory = new File(Properties.reportDirectory, "neighborhoods");
		File indexFile = new File(neighborhoodDirectory, index + Properties.reportWebSuffix);
		
		//find all matches that are from the same chromosome
		ArrayList<SpectrumPeptideMatch> theseMatches = new ArrayList<SpectrumPeptideMatch>();
		for (int i = 0; i < matches.size(); i++) {
			SpectrumPeptideMatch thisMatch = matches.get(i);
			if (thisMatch.getSequence() != match.getSequence()) continue;
			if (Math.abs(thisMatch.getPeptide().getIndex() - match.getPeptide().getIndex()) > Properties.locusNeighborhood) continue;
			theseMatches.add(thisMatch);
		}

		try {
			neighborhoodDirectory.mkdirs();
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(indexFile)));
			appendFile(pw, Properties.reportWebHeaderSubFile);
			
			//print the best match for each spectrum
			pw.println("<h1>Matches in the neighborhood of " + match.getPeptide().getIndex() + "</h1>");
			appendFile(pw, Properties.reportWebTableHeader);
			Collections.sort(theseMatches);
			for (int i = 0; i < theseMatches.size(); i++) {
				pw.println(getTableRow(theseMatches.get(i), i));
			}
			
			appendFile(pw, Properties.reportWebFooterFile);
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
		
		File sequenceDirectory = new File(Properties.reportDirectory, "sequences");
		File indexFile = new File(sequenceDirectory, sequence.getId() + Properties.reportWebSuffix);
		ArrayList<SpectrumPeptideMatch> theseMatches = getMatchesWithSequence(sequence, matches);
		try {
			sequenceDirectory.mkdirs();
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(indexFile)));
			appendFile(pw, Properties.reportWebHeaderSubFile);
			
			//print the best match for each spectrum
			pw.println("<h1>Best match for each spectrum in " + sequence.getSequenceFile().getName() + "</h1>");
			appendFile(pw, Properties.reportWebTableHeader);
			ArrayList<SpectrumPeptideMatch> bestMatches = new ArrayList<SpectrumPeptideMatch>();
			for (int i = 0; i < spectra.size(); i++) {
				ArrayList<SpectrumPeptideMatch> specific = getMatchesWithSpectrum(spectra.get(i), theseMatches);
				Collections.sort(specific);
				if (specific.size() > 0)
					bestMatches.add(specific.get(0));
			}
			Collections.sort(bestMatches);
			for (int i = 0; i < bestMatches.size(); i++) {
				pw.println(getTableRow(bestMatches.get(i), i));
			}
			
			appendFile(pw, Properties.reportWebFooterFile);
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
		File sequenceDirectory = new File(Properties.reportDirectory, "spectra");
		File indexFile = new File(sequenceDirectory, spectrum.getId() + Properties.reportWebSuffix);
		ArrayList<SpectrumPeptideMatch> theseMatches = getMatchesWithSpectrum(spectrum, matches);
		try {
			sequenceDirectory.mkdirs();
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(indexFile)));
			appendFile(pw, Properties.reportWebHeaderSubFile);
			
			//print the best match for each spectrum
			pw.println("<h1>Matches for spectrum " + spectrum.getId() + "</h1>");
			appendFile(pw, Properties.reportWebTableHeader);
			Collections.sort(theseMatches);
			for (int i = 0; i < theseMatches.size(); i++) {
				pw.println(getTableRow(theseMatches.get(i), i));
			}
			
			appendFile(pw, Properties.reportWebFooterFile);
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
	
	private ArrayList<SpectrumPeptideMatch> getMatchesWithSpectrum(Spectrum spectrum, ArrayList<SpectrumPeptideMatch> theseMatches) {
		ArrayList<SpectrumPeptideMatch> out = new ArrayList<SpectrumPeptideMatch>();
		for (int i = 0; i < theseMatches.size(); i++) {
			SpectrumPeptideMatch match = theseMatches.get(i);
			if (match.getSpectrum() == spectrum) {
				out.add(match);
			}
		}
		return out;
	}
	
	private ArrayList<SpectrumPeptideMatch> getMatchesWithSequence(Sequence sequence, ArrayList<SpectrumPeptideMatch> theseMatches) {
		ArrayList<SpectrumPeptideMatch> out = new ArrayList<SpectrumPeptideMatch>();
		for (int i = 0; i < theseMatches.size(); i++) {
			SpectrumPeptideMatch match = theseMatches.get(i);
			if (match.getSequence() == sequence) {
				out.add(match);
			}
		}
		return out;
	}

	private String getTableRow(SpectrumPeptideMatch match, int index) {
		StringBuffer sb = new StringBuffer();
		sb.append("<tr>");
		
		sb.append("<td>");
		sb.append("<a href=\"../spectra/");
		sb.append(match.getSpectrum().getId());
		sb.append(Properties.reportWebSuffix);
		sb.append("\">");
		sb.append(match.getSpectrum().getId());
		sb.append("</a>");
		sb.append("</td>");
		
		sb.append("<td>");
		sb.append(match.getPeptide().getSequence());
		sb.append("</td>");
		
		sb.append("<td><nobr>");
		sb.append("<a href=\"../sequences/");
		sb.append(match.getSequence().getId());
		sb.append(Properties.reportWebSuffix);
		sb.append("\">");
		sb.append(match.getSequence().getSequenceFile().getName ());
		sb.append("</a>");
		if (match.getPeptide().isForward()) {sb.append(" forward ");}
		else {sb.append(" reverse ");}
		sb.append(match.getPeptide().getReadingFrame());
		sb.append("</nobr></td>");
		
		sb.append("<td>");
		sb.append(match.getPeptide().getIndex());
		sb.append("</td>");
		
		sb.append("<td>");
		sb.append(match.getScoreHMM());
		sb.append("</td>");
		
		sb.append("<td>");
		sb.append(match.getMSMSFitRank());
		sb.append("</td>");
		
		return sb.toString();
	}
	

	private void appendFile(PrintWriter pw, File file) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			String line = br.readLine();
			while (line != null) {
				pw.println(line);
				line = br.readLine();
			}
			br.close();
		} catch (FileNotFoundException e) {
			U.p("could not append file: " + file.getName());
			e.printStackTrace();
		} catch (IOException e) {
			U.p("could not read file: " + file.getName());
			e.printStackTrace();
		}
	}

}
