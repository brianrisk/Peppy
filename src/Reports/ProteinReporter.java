package Reports;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import Peppy.Properties;
import Peppy.Protein;
import Utilities.U;

public class ProteinReporter {
	
	private ArrayList<Protein> proteins;
	
	private PrintWriter pw;
	
//	public static void main(String args[]) {
//		print("<td bgcolor=\"#" + getHexColor(Color.red) + "\">");
//	}
	
	public ProteinReporter(ArrayList<Protein> proteins, File destinationFile) {
		this.proteins = proteins;
		try {
			pw = new PrintWriter(new BufferedWriter(new FileWriter(destinationFile)));
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void generateReport() {
		U.appendFile(pw, Properties.reportWebHeaderFile);
		for (Protein protein: proteins) {
			printProtein(protein, pw);
		}
		U.appendFile(pw, Properties.reportWebFooterFile);
		pw.flush();
		pw.close();
	}
	
	private void printProtein(Protein protein, PrintWriter pw) {
		String acidString = protein.getAcidString();
		int [] matchPositions = protein.getMatchPositions();
		printH3(protein.getName());
		printP("Match coverage: " + protein.getMatchCoverage());
		printP("Protein score: " + protein.getScore());
		print("<table>");
		
		Color PINK = new Color(255, 128, 128);
		for (int i = 0; i < acidString.length(); i++) {
			if (i % 20 == 0) {
				print("<tr>");
				printTD("" + i);
			}
			if (matchPositions[i] == Protein.T_FPR01) printTD("" + acidString.charAt(i), Color.RED);
			if (matchPositions[i] == Protein.T_FPR05) printTD("" + acidString.charAt(i), PINK);
			if (matchPositions[i] == Protein.T_FPRXX) printTD("" + acidString.charAt(i), Color.LIGHT_GRAY);
			if (matchPositions[i] == Protein.T_MOD) printTD("" + acidString.charAt(i), Color.GREEN);
			if (matchPositions[i] == Protein.T_NOTHING) printTD("" + acidString.charAt(i));
		}
		print("</table>");
		
	}
	
	private void printP(String string) {
		print("<p>");
		print(string);
		print("</p>");
	}
	
	private void printH1(String string) {printH(string, 1);}
	private void printH2(String string) {printH(string, 2);}
	private void printH3(String string) {printH(string, 3);}
	
	private void printH(String string, int level) {
		print("<h" + level + ">");
		print(string);
		print("</h" + level + ">");
	}
	
	private void printTD(String string) {
		print("<td>");
		print(string);
		print("</td>");
	}
	private void printTD(String string, Color color) {
		print("<td bgcolor=\"#" + getHexColor(color) + "\">");
		print(string);
		print("</td>");
	}
	
	private void print(String string) {
		pw.println(string);
	}
	
	private static String getHexColor(Color color) {
		String out = Integer.toHexString(color.getRGB());
		out = out.substring(2, out.length());
		return out;
	}

}
