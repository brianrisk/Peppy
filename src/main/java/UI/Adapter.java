package UI;

import processing.core.PApplet;

import java.util.ArrayList;

public abstract class Adapter {

    ArrayList<View> views = new ArrayList<View>();

    public void mousePressed(PApplet context) {
        for (int index = 0; index < views.size(); index++) {
            View theView = (View) views.get(index);
            if (theView.isActive) {
                if (theView.isIn(context.mouseX, context.mouseY)) theView.mousePressed();
            }
        }
    }

    public void mouseDragged(PApplet context) {
        for (int index = 0; index < views.size(); index++) {
            View theView = (View) views.get(index);
            if (theView.isActive) {
                if (theView.isIn(context.mouseX, context.mouseY)) theView.mouseDragged();
            }
        }
    }

    public void mouseReleased(PApplet context) {
        for (int index = 0; index < views.size(); index++) {
            View theView = (View) views.get(index);
            if (theView.isActive) {
                if (theView.isIn(context.mouseX, context.mouseY)) theView.mouseReleased();
            }
        }
    }

    public void render() {
        for (int index = 0; index < views.size(); index++) {
            View theView = (View) views.get(index);
            if (theView.isActive) {
                theView.render();
            }
        }
    }
}