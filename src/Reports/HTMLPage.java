package Reports;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import Peppy.Properties;
import Utilities.U;

public abstract class HTMLPage {
	

	protected File destinationFile;
	private PrintWriter pw;
	
	
	public HTMLPage(File destinationFile) {
		this.destinationFile = destinationFile;
		try {
			pw = new PrintWriter(new BufferedWriter(new FileWriter(destinationFile)));
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public abstract void makePage();
	
	
	protected void printP(String string) {
		print("<p>");
		print(string);
		print("</p>");
	}
	
	protected void printH1(String string) {printH(string, 1);}
	
	protected void printH2(String string) {printH(string, 2);}
	
	protected void printH3(String string) {printH(string, 3);}
	
	protected void printH(String string, int level) {
		print("<h" + level + ">");
		print(string);
		print("</h" + level + ">");
	}
	
	protected void printTD(String string) {
		print("<td>");
		print(string);
		print("</td>");
	}
	protected void printTD(String string, Color color) {
		print("<td bgcolor=\"#" + getHexColor(color) + "\">");
		print(string);
		print("</td>");
	}
	
	protected void printTH(String string) {
		print("<th>");
		print(string);
		print("</th>");
	}
	
	protected void printLI(String string) {
		print("<li>");
		print(string);
		print("</li>");
	}
	
	protected void print(String string) {
		pw.println(string);
	}
	
	protected void printHeader() {
		U.appendFile(pw, Properties.reportWebHeaderFile);
	}
	
	protected void printFooter() {
		U.appendFile(pw, Properties.reportWebFooterFile);
		pw.flush();
		pw.close();
	}
	
	protected static String getHexColor(Color color) {
		String out = Integer.toHexString(color.getRGB());
		out = out.substring(2, out.length());
		return out;
	}
	

}
