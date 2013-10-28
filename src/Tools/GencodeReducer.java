package Tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import Peppy.U;

public class GencodeReducer {
	
	public static void main(String args[]) {
		
		
		try {
			File gencodeFile = new File("resources/gencode/gencode.v16.annotation.gtf");
			BufferedReader gencodeReader = new BufferedReader(new FileReader(gencodeFile));
			
			File gencodeReducedFile = new File("gencodeReduced.gtf");
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(gencodeReducedFile)));
			
			String line = gencodeReader.readLine();
			while (line != null) {
				if (line.indexOf("\ttranscript\t") != -1) {
//					if (line.indexOf("transcript_type \"protein_coding\"") != -1) {
						pw.println(line);
//					}
				}
				line = gencodeReader.readLine();
			}
			
			gencodeReader.close();
			pw.flush();
			pw.close();
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		U.p("done");
	}

}
