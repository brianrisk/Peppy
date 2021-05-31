package UI;

import Peppy.U;
import processing.core.PApplet;
import processing.core.PFont;

import java.lang.management.ManagementFactory;
import java.lang.management.OperatingSystemMXBean;

public class UserInterface extends PApplet {

	private static final long serialVersionUID = 5616994580457892732L;
	private Stack stack;

	public static void main(String args[]) {
		PApplet.main(new String[] {"UI.UserInterface" });	
	}

	public void settings() {
		size(1024, 768, "processing.awt.PGraphicsJava2D");
		// TODO: commenting pixelDensity for now as it clips the menu items
//		pixelDensity(2);
	}

	public void setup() {
//		smooth(8);
		UIC.init(this);
		UIF.init(this);
		frameRate(200);
		textFont(UIC.arial, UIC.fontSize);
		
		stack = new Stack(this);
	}

	public void draw() {
		stack.render();
		OperatingSystemMXBean operatingSystemMXBean = ManagementFactory.getOperatingSystemMXBean();
		double load = operatingSystemMXBean.getSystemLoadAverage() / operatingSystemMXBean.getAvailableProcessors();
		fill(255,255,255,255);
		UIF.drawText(UIC.gridX * 14, UIC.gridY * 0, UIC.gridX * 2, "CPU: " + load, UIC.ALIGN_RIGHT);
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
