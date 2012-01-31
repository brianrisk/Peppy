package Tools;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import Peppy.U;


public class CleanCode {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		File codeFolder = new File("JainabCode");
		File cleanedFolder = new File("JainabCleaned");
		File[] files = codeFolder.listFiles();
		for (int i = 0; i < files.length; i++) {
			if (!files[i].getName().endsWith(".java")) continue;
			if(files[i].isHidden()) continue;
			try {
				BufferedReader br = new BufferedReader(new FileReader(files[i]));
				PrintWriter pw = new PrintWriter(new FileWriter(new File(cleanedFolder, files[i].getName())));
				String line = br.readLine();
				while (line != null) {
					line = line.trim();
					if (line.startsWith("//")) {line = br.readLine(); continue;}
					if (line.startsWith("System.out")) {line = br.readLine(); continue;}
					pw.println(line);
					line = br.readLine();
				}
				pw.flush();
				pw.close();
				br.close();
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		U.p("done a little spring cleaning.");

	}

}
