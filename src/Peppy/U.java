package Peppy;

import java.awt.Color;
import java.awt.Toolkit;
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
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.Stack;

/**
 * A basic utility class.
 * Saves a lot of time not writing out "System.out.println"
 * 
 * Copyright 2012, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class U {
	public static final long SECOND = 1000;
	public static final long MINUTE = SECOND * 60;
	public static final long HOUR = MINUTE * 60;
	public static final long DAY = HOUR * 24;
	public static final long YEAR = DAY * 365;
	

	private static Stack<Long> stopwatchClicks = new Stack<Long>();
	
	public static void startStopwatch() {
		stopwatchClicks.push(System.currentTimeMillis());
	}
	
	public static long stopStopwatch() {
		long stopTimeMilliseconds = System.currentTimeMillis();
		long timeElapsed = stopTimeMilliseconds - stopwatchClicks.pop();
		U.p("Time Elapsed: " + U.millisecondsToString(timeElapsed));
		return timeElapsed;
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
	
	public static String getDateString() {
		String response = "";
		GregorianCalendar calendar = new GregorianCalendar();
		Date time = new Date();
		calendar.setTime(time);
		response += calendar.get(GregorianCalendar.YEAR);
		response += "-";
		response += calendar.get(GregorianCalendar.MONTH);
		response += "-";
		response += calendar.get(GregorianCalendar.DAY_OF_MONTH);
		response += " ";
		response += calendar.get(GregorianCalendar.HOUR_OF_DAY);
		response += "-";
		if (calendar.get(GregorianCalendar.MINUTE) < 10) {
			response += "0";
		}
		response += calendar.get(GregorianCalendar.MINUTE);
		return response;
	}
	
	public static void p(Object o) {System.out.println(o);}
	public static void p(double o) {System.out.println(o);}
	public static void p(int o) {System.out.println(o);}
	public static void p(char o) {System.out.println(o);}
	public static void p() {System.out.println();}
	
	public static String in() {
		BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
		String out = "";
		try {
			out = br.readLine();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return out;
	}
	
	public static void beep() {Toolkit.getDefaultToolkit().beep();}
	
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
	
	public static double log(double base, double of) {
		return Math.log(of) / Math.log(base);
	}
	
	public static String getFileNameWithoutSuffix(File file) {
        return file.getName().substring(0, file.getName().lastIndexOf('.'));
	}
	
	public static String reverseString(String in) {
		return new StringBuffer(in).reverse().toString();
	}
	
	public static String getRGBStringFromPercent(double percent) {
		Color hsb = Color.getHSBColor((float) percent, 1.0f, 1.0f);
		String rgb = Integer.toHexString(hsb.getRGB());
		rgb = rgb.substring(2, rgb.length());
		return rgb;
	}

}
