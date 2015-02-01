package UI;

import processing.core.PApplet;

public class UserInterface extends PApplet {

	private static final long serialVersionUID = 5616994580457892732L;

	public static void main(String args[]) {
		PApplet.main(new String[] {"UI.UserInterface" });
	}

	public void setup() {
		size(200,200);
		background(0);
	}

	public void draw() {
		stroke(255);
		if (mousePressed) {
			line(mouseX,mouseY,pmouseX,pmouseY);
		}
	}
}
