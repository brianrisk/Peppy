package Reports;

/**
 * generates some useful header, table, footer strings and the like
 * @author Brian Risk
 *
 */
public class ReportStrings {
	
	public static String getHeader() {
		StringBuffer sb = new StringBuffer();
		sb.append("<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">");
		sb.append("<html xmlns=\"http://www.w3.org/1999/xhtml\">");
		sb.append("<head>");
		sb.append("<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\" />");
		sb.append("<title>Peppy Report</title>");
		sb.append("<style type=\"text/css\">");
		sb.append("<!--");
		sb.append("@import url(\"http://proteomics.me/resources/reports/style.css\");");
		sb.append("@import url(\"http://proteomics.me/resources/reports/sortable.css\");");
		sb.append("-->");
		sb.append("</style>");
		sb.append("<script type=\"text/javascript\" src=\"http://proteomics.me/resources/reports/sortable.js\"></script>");
		sb.append("</head>");
		return sb.toString();
	}
	
	public static String getTableHeader() {
		StringBuffer sb = new StringBuffer();
		sb.append("<table class=\"sortable\" id=\"box-table-a\">");
		sb.append("    <thead>");
		sb.append("    	<tr>");
		sb.append("            <th scope=\"col\">ID</th>");
		sb.append("            <th scope=\"col\">Peptide</th>");
		sb.append("            <th scope=\"col\">Chromosome</th>");
		sb.append("            <th scope=\"col\">Start</th>");
		sb.append("            <th scope=\"col\">Stop</th>");
		sb.append("            <th scope=\"col\">F</th>");
		sb.append("            <th scope=\"col\">Spliced</th>");
		sb.append("            <th scope=\"col\">Mods</th>");
		sb.append("            <th scope=\"col\">Score</th>");
		sb.append("            <th scope=\"col\">ions</th>");
		sb.append("            <th scope=\"col\">ions/length</th>");
		sb.append("            <th scope=\"col\">score ratio</th>");
		sb.append("            <th scope=\"col\">E value</th>");
		sb.append("        </tr>");
		sb.append("    </thead>");
		sb.append("     <tbody>");
		return sb.toString();
	}
	
	public static String getFooter() {
		StringBuffer sb = new StringBuffer();
		sb.append("</tbody>");
		sb.append("</table>");
		sb.append("<body>");
		sb.append("<html>");
		return sb.toString();
	}
	

}
