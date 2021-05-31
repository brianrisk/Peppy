package Tools;

import Peppy.U;

import java.io.*;

public class GencodeReducer {
	
	public static void main(String args[]) {
		
		
		try {
			File gencodeFile = new File("resources/gencode/gencode.v16.annotation.gtf");
			BufferedReader gencodeReader = new BufferedReader(new FileReader(gencodeFile));
			
			File gencodeReducedFile = new File("gencodeUTRReduced.gtf");
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(gencodeReducedFile)));
			
			String line = gencodeReader.readLine();
			while (line != null) {
				
				//get the full transcript length
//				if (line.indexOf("\ttranscript\t") != -1) {
				if (line.indexOf("\tUTR\t") != -1) {
					
					//use only protein-coding transcripts
					if (line.indexOf("transcript_type \"protein_coding\"") != -1) {
						
						pw.println(line);
					}
				}
				line = gencodeReader.readLine();
			}
			
			gencodeReader.close();
			pw.flush();
			pw.close();
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		U.p("done");
	}

}
