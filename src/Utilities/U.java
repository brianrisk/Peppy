package Utilities;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintWriter;

/**
 * A basic utility class.
 * Saves a lot of time not writing out "System.out.println"
 * @author Brian Risk
 *
 */
public class U {
	public static final long SECOND = 1000;
	public static final long MINUTE = SECOND * 60;
	public static final long HOUR = MINUTE * 60;
	public static final long DAY = HOUR * 24;
	public static final long YEAR = (long) (DAY * 365.25);
	
	private static long startTimeMilliseconds;
	private static long stopTimeMilliseconds;
	
	public static void startStopwatch() {
		startTimeMilliseconds = System.currentTimeMillis();
	}
	
	public static long stopStopwatch() {
		stopTimeMilliseconds = System.currentTimeMillis();
		long timeElapsed = stopTimeMilliseconds - startTimeMilliseconds;
		U.p("Time Elapsed: " + U.millisecondsToString(timeElapsed));
		return timeElapsed;
	}
	
	public static void printTimeRemaining(double amountComplete) {
		double elapsed = System.currentTimeMillis() - startTimeMilliseconds;
		long timeRemaining = (long) ((elapsed / amountComplete) - elapsed);
		U.p("Time Remaining: " + U.millisecondsToString(timeRemaining));
	}
	
	public static void save(String fileName, String output) {
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(fileName)));
			pw.print(output);
			pw.flush();
			pw.close();
		}catch (IOException ioe) {
			p("could not save file: " + fileName);
			ioe.printStackTrace();
		}
	}
	
	public static String readFileToString(String fileName) {
		StringBuffer sb = new StringBuffer();
		try {
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			String line = br.readLine();
			while (line != null) {
				sb.append(line);
				line = br.readLine();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return sb.toString();
	}
	
	public static void printUserDirectory() {
		String tmp=System.getProperty("user.dir");
	    U.p(tmp);
	}
	
	public static String millisecondsToString(long timeElapsed) {
		String response = "";
		long amount;
		if (timeElapsed >= YEAR) {
			amount = timeElapsed / YEAR;
			response += amount + " year";
			if (amount != 1) response += "s";
			response += ", ";
			timeElapsed -= YEAR * amount;
		}
		if (timeElapsed >= DAY) {
			amount = timeElapsed / DAY;
			response += amount + " day";
			if (amount != 1) response += "s";
			response += ", ";
			timeElapsed -= DAY * amount;
		}
		if (timeElapsed >= HOUR) {
			amount = timeElapsed / HOUR;
			response += amount + " hour";
			if (amount != 1) response += "s";
			response += ", ";
			timeElapsed -= HOUR * amount;
		}
		if (timeElapsed >= MINUTE) {
			amount = timeElapsed / MINUTE;
			response += amount + " minute";
			if (amount != 1) response += "s";
			response += ", ";
			timeElapsed -= MINUTE * amount;
		}
		if (timeElapsed >= SECOND) {
			amount = timeElapsed / SECOND;
			response += amount + " second";
			if (amount != 1) response += "s";
			timeElapsed -= SECOND * amount;
		}
		return response;
	}
	
	public static void p(Object o) {System.out.println(o);}
	public static void p(double o) {System.out.println(o);}
	public static void p(int o) {System.out.println(o);}
	public static void p(char o) {System.out.println(o);}
	public static void p() {System.out.println();}
	
	public static void copyfile(File sourceFile, File destinationFile){
	    try{
	      InputStream in = new FileInputStream(sourceFile);
	      OutputStream out = new FileOutputStream(destinationFile);

	      byte[] buffer = new byte[1024];
	      int len;
	      while ((len = in.read(buffer)) > 0){
	        out.write(buffer, 0, len);
	      }
	      in.close();
	      out.close();
	    }
	    catch(FileNotFoundException ex){
	      System.out.println(ex.getMessage() + " in the specified directory.");
	      System.exit(0);
	    }
	    catch(IOException e){
	      System.out.println(e.getMessage());      
	    }
	  }
	
	
	public static void appendFile(PrintWriter pw, File file) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			String line = br.readLine();
			while (line != null) {
				pw.println(line);
				line = br.readLine();
			}
			br.close();
		} catch (FileNotFoundException e) {
			U.p("could not append file: " + file.getName());
			e.printStackTrace();
		} catch (IOException e) {
			U.p("could not read file: " + file.getName());
			e.printStackTrace();
		}
	}
	
	public static double log(double base, double of) {
		return Math.log(of) / Math.log(base);
	}
	
	public static double calculateM(double [] xValues, double [] yValues, int start, int stop) {
		double numerator1, numerator2, denomenator1, denomenator2;
		double numerator = 0.0, denomenator = 0.0;
		double temp1 = 0.0, temp2 = 0.0, temp = 0.0;
		double parameterM;
		int i;
		for (i = start; i < stop; i++) {
			temp += (xValues[i] * yValues[i]);
		}
		numerator1 = (stop - start) * (temp);
		for (i = start; i < stop; i++) {
			temp1 += xValues[i];
		}
		for (i = start; i < stop; i++) {
			temp2 += yValues[i];
		}
		numerator2 = temp1 * temp2;
		numerator = numerator1 - numerator2;
		temp1 = temp2 = 0.0;
		for (i = start; i < stop; i++) {
			temp1 = xValues[i];
			temp2 += (temp1 * temp1);
		}
		denomenator1 = (stop - start) * temp2;

		temp1 = 0.0; 
		for (i = start; i < stop; i++) {
			temp1 += xValues[i];
		}
		denomenator2 = (temp1 * temp1);
		denomenator = denomenator1 - denomenator2;
		parameterM = numerator / denomenator;
		return parameterM;
	}
	
	public static double calculateB(double [] xValues, double [] yValues, int start, int stop, double m) {
		double parameterB;
		double temp1, temp2;
		int i;

		temp1 = temp2 = 0.0;
		for (i = start; i < stop; i++) {
			temp1 += xValues[i];
		}

		for (i = start; i < stop; i++) {
			temp2 += yValues[i];
		}
		parameterB = (1.0 / (stop - start)) * (temp2 - m * temp1);
		return parameterB;
	}

}
