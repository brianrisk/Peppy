package Reports;

import java.io.File;

public class Properties {
	
	//Report related
	public static File reportDirectory = new File("reports");
	public static String reportWebSuffix = ".html";
	public static File reportWebHeaderFile = new File("resources/reports/header.txt");
	public static File reportWebHeaderSubFile = new File("resources/reports/header-sub.txt");
	public static File reportWebFooterFile = new File("resources/reports/footer.txt");
	public static File reportWebTableHeader = new File("resources/reports/index-table-header.txt");
	
	//the number of nucleotides away from a specific location on a chromosome for it to be
	//considered part of the "neighborhood"
	public static int locusNeighborhood = 3000;

}
