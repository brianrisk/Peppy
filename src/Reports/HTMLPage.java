package Reports;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Copyright 2012, Brian Risk
 * Released under the Netscape Public License
 * 
 * @author Brian Risk
 *
 */
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
	
	protected void printTR() {
		print("<tr>");
	}
	
	protected void printLI(String string) {
		print("<li>");
		print(string);
		print("</li>");
	}
	
	protected void printLink(String link, String text) {
		print("<a href=\"" + link + "\">" + text + "</a>");
	}
	
	protected void print(String string) {
		pw.println(string);
	}
	
	protected void printHeader() {
		printHeader("Peppy Report", "");
	}
	
	protected void printHeader(String title, String other) {
		print("<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">");
		print("<html xmlns=\"http://www.w3.org/1999/xhtml\">");
		print("<head>");
		print("<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\" />");
		print("<title>" + title + "</title>");
		print(other);
		print("<style type=\"text/css\">");
		print("<!--");
		print("@import url(\"http://proteomics.me/resources/reports/style.css\");");
		print("@import url(\"http://proteomics.me/resources/reports/sortable.css\");");
//		print("@import url(\"../../style.css\");");
//		print("@import url(\"../../sortable.css\");");
		print("-->");
		print("</style>");
		print("<script type=\"text/javascript\" src=\"http://proteomics.me/resources/reports/sortable.js\"></script>");
//		print("<script type=\"text/javascript\" src=\"../../sortable.js\"></script>");
		print("</head>");
		//title
		print("<div id=\"container\">");
		print("<div id=\"header\">");
		printH1(title);
		print("</div>");
		print("<div id=\"emphasis-bar-top\"><p></p></div>");
		print("<div id=\"body\">");
	}
	
	protected void printFooter() {
		print("</div>");
		print("<div id=\"footer\">");
		print("<div id=\"emphasis-bar\"><p></p></div>");
		print("<p>Created with <a href=\"http://peppyresearch.com/peppy\">Peppy</a>, protein identification, proteogenomic mapping software ");
//		print("by <a href=\"http://unitedproteomics.com\">United Proteomics</a></p>");
		print("</div>");
		print("</div>");
		print("</body>");
		print("</html>");
		pw.flush();
		pw.close();
	}
	
	protected static String getHexColor(Color color) {
		String out = Integer.toHexString(color.getRGB());
		out = out.substring(2, out.length());
		return out;
	}
	

}
