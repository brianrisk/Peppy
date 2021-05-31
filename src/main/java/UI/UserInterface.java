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
//		hint(ENABLE_RETINA_PIXELS);
//		smooth(8);
		size(1024,768,JAVA2D);
		
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
