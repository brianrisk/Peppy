package UI;

import java.lang.management.ManagementFactory;
import java.lang.management.OperatingSystemMXBean;

import processing.core.PApplet;
import processing.core.PFont;
import Peppy.U;

public class UserInterface extends PApplet {

	private static final long serialVersionUID = 5616994580457892732L;
	private Stack stack;

	public static void main(String args[]) {
		PApplet.main(new String[] {"UI.UserInterface" });	
	}

	public void setup() {
		size(1024,768,PApplet.JAVA2D);
		Constants.init(this);
		Functions.init(this);
		textFont(Constants.arial, Constants.fontSize);
		smooth();
		frameRate(200);
		stack = new Stack(this);
	}

	public void draw() {
		stack.render();
		OperatingSystemMXBean operatingSystemMXBean = ManagementFactory.getOperatingSystemMXBean();
		double load = operatingSystemMXBean.getSystemLoadAverage() / operatingSystemMXBean.getAvailableProcessors();
		fill(255,255,255,255);
		Functions.drawText(Constants.gridX * 14, Constants.gridY * 0, Constants.gridX * 2, "FPS: " + load, Constants.ALIGN_RIGHT);
	}
	

	public void mousePressed() {
		stack.mousePressed(this);
	}

	public void mouseDragged() {
		stack.mouseDragged(this);
	}

	public void mouseReleased() {
		stack.mouseReleased(this);
	}
	
	public static void printFonts() {
		String [] fonts = PFont.list();
		for (String font: fonts) {
			U.p(font);
		}
	}
	

}
