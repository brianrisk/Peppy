package Experimental;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import Peppy.U;

public class ScottNicholson {
	
	//scott's scaffold file came as one big file.
	//this seems to be causing some problems.
	//Let's break it up.
	public static void main(String [] args) {
		File genome = new File("/Volumes/Research/scott-nicholson/sequences/pea_aphid_genome.1.fasta");
		File saveDir = new File ("/Volumes/Research/scott-nicholson/sequences/scaffolds");
		saveDir.mkdirs();
		try {
			int contigIndex = 0;
			BufferedReader reader = new BufferedReader(new FileReader(genome));
			String line = reader.readLine();
			PrintWriter writer = null;
			while (line != null) {
				if (line.startsWith(">")) {
					if (writer != null) {
						writer.flush();
						writer.close();
					}
					contigIndex++;
					writer = new PrintWriter(new BufferedWriter(new FileWriter(new File(saveDir, "scaffold-" + contigIndex + ".fasta"))));
				}
				writer.println(line);
				line = reader.readLine();
			}
		} catch (FileNotFoundException e) {

			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		U.p("done");
	}

}
