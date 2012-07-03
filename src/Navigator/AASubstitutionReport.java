package Navigator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;

import Math.MassError;
import Peppy.Properties;
import Peppy.U;

/**
 * Copyright 2012, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class AASubstitutionReport {
	
	public static void main(String args[]) {
		new AASubstitutionReport();
		U.p("done");
	}
	
	public AASubstitutionReport() {
		ArrayList<AASubstitution> substitutions = AASubstitution.generateListOfAASubstitutions();
		
		/* score cutoff for loads samples */
		double matchScoreCutoff = 14.24;
		
		/* load the sample */
		Sample whim2Human = new Sample("21881_FA_Example1_dta", Sample.REFERENCE_PROTEIN);
//		whim2Human.loadResults(new File("/Users/risk2/PeppyData/Mayo/reports/21881_FA_Example1_dta/report.txt"));
		whim2Human.loadResults(new File("/Users/risk2/PeppyData/Mayo/reports/Mayo dogan protein/report.txt"));
		
		whim2Human.reduceMatches(matchScoreCutoff);
		
		try {
			PrintWriter matchWriter;
			matchWriter = new PrintWriter(new BufferedWriter(new FileWriter(new File("AASubstitutions.html"))));
		
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
		
			ArrayList<Match> matches = whim2Human.getMatches();
			for (Match match: matches) {
				if (match.getBoolean("isModified")) {
					String sequenceName = match.getString("sequenceName");
					double modMass = match.getDouble("modMass");
					int modIndex = match.getInt("modIndex");
					double peptideMass = match.getDouble("PrecursorNeutralMass");
					String peptide = match.getString("peptideSequence");
					char AA = peptide.charAt(modIndex);
					
					/* ignore oxidations */
					if (AA == 'M' && Math.round(modMass * 10) == 160) continue;
						
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

}
