package Experimental;

import Peppy.U;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class ProteoWizardPathMaker {
	
	public static void main(String args []) {
		File source = new File ("/Volumes/Research/Breast-mzml");
		File destination = new File ("/Volumes/Research/Breast-converted");
		
		try {
			PrintWriter pw = new PrintWriter(new FileWriter(new File("/Volumes/Research/workspace/Peppy/proteowizard.sh")));
			root(source, destination, "", pw);
			pw.flush();
			pw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	/* recursively go into diretories */
	public static void root(File source, File destination, String path, PrintWriter pw) {
		
		if (source.isDirectory()) {
			File [] sourceFiles = source.listFiles();
			for (File sourceFile: sourceFiles) {
				root(sourceFile, destination, path + source.getName() + "/", pw);
			}
		}
		else {
			if (source.getName().toLowerCase().endsWith(".mzml")) {
				File saveDir = new File(destination, path);
				saveDir.mkdirs();
				U.p("./msconvert -o " + saveDir.getAbsolutePath() + " " + source.getAbsolutePath() + " --mgf");
				pw.println("./msconvert -o " + saveDir.getAbsolutePath() + " " + source.getAbsolutePath() + " --mgf");
			}
		}
		
	}

}
