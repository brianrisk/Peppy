package Navigator;

import Math.MassError;
import Peppy.Properties;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Copyright 2013, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class AASubstitutionReport {
	
	File variModFile;
	ArrayList<Match> matches;
	

	/**
	 * 
	 * Example:
	 * new AASubstitutionReport(new File ("varimod/report.txt"), new File("AASubstitutions.html"));
	 * 
	 * @param variModFile The .txt report file for the variable modifications
	 * @param saveFile .html file where this will be saved
	 */
	public AASubstitutionReport(File variModFile, File saveFile) {
		this.variModFile = variModFile;
		ArrayList<AASubstitution> substitutions = AASubstitution.generateListOfAASubstitutions();
		
		/* load the matches */
		matches = Match.loadMatches(variModFile);	
		
		try {
			PrintWriter matchWriter;
			matchWriter = new PrintWriter(new BufferedWriter(new FileWriter(saveFile)));
		
			matchWriter.println(HTML.sortableTableHeader);
			matchWriter.println(HTML.tableTop);
			matchWriter.println("<tr>");
			matchWriter.println("<th>sequenceName</th>");
			matchWriter.println("<th>modIndex</th>");
			matchWriter.println("<th>previous</th>");
			matchWriter.println("<th>present</th>");
			matchWriter.println("<th>peptide</th>");
			matchWriter.println("<th>modMass</th>");
			matchWriter.println("</tr>");
			matchWriter.println("</thead>");
			matchWriter.println("<tbody>");
		
			for (Match match: matches) {
				if (match.getBoolean("isModified") && match.getBoolean("modLocCertain")) {
					String sequenceName = match.getString("sequenceName");
					double modMass = match.getDouble("modMass");
					int modIndex = match.getInt("modIndex");
					double peptideMass = match.getDouble("PrecursorNeutralMass");
					String peptide = match.getString("peptideSequence");
					char AA = peptide.charAt(modIndex);
					
					/* trim the sequence name if it has a pipe character */
					if (sequenceName.indexOf("|") != -1) {
						sequenceName = sequenceName.split("\\|")[0];
					}
					
					/* ignore oxidation */
//					if (AA == 'M' && Math.round(modMass) == 16) continue;
					if (Math.round(modMass) == 16) continue;
					
					/* ignore deamidation */
					if (AA == 'N' && Math.round(modMass) == 1) continue;
					
					/* ignore carbamylation */
//					if ((AA == 'K'  || AA == 'R'  || AA == 'C'  || modIndex == 0) && Math.round(modMass) == 43) continue;
					if (Math.round(modMass) == 43) continue;
						
					double error = MassError.getDaltonError(Properties.fragmentTolerance, peptideMass);
					int start = Math.abs(Collections.binarySearch(substitutions, new AASubstitution(modMass - error))) - 1;
					int stop = Math.abs(Collections.binarySearch(substitutions, new AASubstitution(modMass + error))) - 1;
					for (int i = start; i < stop; i++) {
						if (substitutions.get(i).getPrevious() == AA) {
							
							
							
							matchWriter.println("<tr>");
							matchWriter.println("<td>" + sequenceName + "</td>");
							matchWriter.println("<td>" + modIndex + "</td>");
							matchWriter.println("<td>" + AA + "</td>");
							matchWriter.println("<td>" + substitutions.get(i).getPresent() + "</td>");
							matchWriter.println("<td>" + peptide + "</td>");
							matchWriter.println("<td>" + modMass + "</td>");
							matchWriter.println("</tr>");
							break;
						}
					}		
				}
				
			}
			matchWriter.println(HTML.sortableTableFooter);
			matchWriter.flush();
			matchWriter.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public ArrayList<Match> getMatches() {
		return matches;
	}

}
