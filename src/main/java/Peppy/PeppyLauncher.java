package Peppy;

import java.io.IOException;
import java.net.URISyntaxException;

/**
 * Launches the main Peppy class and allocates the maximum amount of physical memory possible
 * 
 * References:
 * http://silentdevelopment.blogspot.com/2010/03/how-to-set-or-increase-xmx-heap-memory.html
 * http://alvinalexander.com/java/java-exec-processbuilder-process-1
 * 
 * @author brianrisk
 *
 */
public class PeppyLauncher {
	
	public static void main(String [] args) {
		try {
			String pathToJar = PeppyLauncher.class.getProtectionDomain().getCodeSource().getLocation().toURI().getPath();
			U.p(pathToJar);
			ProcessBuilder pb = new ProcessBuilder("java","-XX:+AggressiveHeap", "-classpath", pathToJar, "Peppy.Peppy");
			pb.start();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (URISyntaxException e) {
			e.printStackTrace();
		}
	}

}
