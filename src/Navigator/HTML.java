package Navigator;

/**
 * 
 * Static methods returning HTML strings
 * 
 * Copyright 2013, Brian Risk
 * 
 * @author Brian Risk
 *
 */
public class HTML {
	
	public static String sortableTableHeader = 
					"<html>" +
					"<head>" +
					"<style type=\"text/css\">" +
					"<!--" +
					"@import url(\"http://geneffects.com/files/peppy/css/blue-table/style.css\");" +
					"-->" +
					"</style>" +
					"<script type=\"text/javascript\" src=\"http://geneffects.com/files/peppy/js/jquery-1.7.2.min.js\"></script>" +
					"<script type=\"text/javascript\" src=\"http://geneffects.com/files/peppy/js/jquery.tablesorter.min.js\"></script>" +
					"<script type=\"text/javascript\">" +
					"$(document).ready(function() " +
					"    { " +
					"        $(\"#myTable\").tablesorter(); " +
					"    } " +
					"); " +
					"</script>" +
					"</head>" +
					"<body>";
	
	public static String sortableTableFooter = "</tbody></table></body></html>";
	
	public static String tableTop = "<table class=\"tablesorter\" id=\"myTable\"><thead>";

}
