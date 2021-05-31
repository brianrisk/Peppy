package UI;

import processing.core.PApplet;
import processing.core.PGraphics;
import processing.core.PImage;

public class ViewList extends View {

    PApplet context;

    String[] labels;

    int selectedIndex = 0;
    boolean hasShadow = true;

    PGraphics slide;
    int slideWidth;
    int slideHeight;
    int slideY = 0;
    int slideYMax;
    int slideYMin;
    int slideYOld = 0;
    int mousePressedY = 0;
    int shadeHeight;

    /************************************
     * Methods to override / implement
     ************************************/

    public void render() {
        renderList();
    }


    public void updateLocation(int x, int y) {
        viewX = x;
        viewY = y;
    }


    /************************************
     * Public methods
     ************************************/


    /*
     * Constructor
     */
    public ViewList(PApplet context, Listener listener, int x, int y, int width, int height) {
        this.context = context;
        this.listener = listener;
        this.viewX = x;
        this.viewY = y;
        this.viewWidth = width;
        this.viewHeight = height;

        //the shading
        shadeHeight = 1 * UIC.gridY;

    }

    public void setLabels(String[] labels) {
        this.labels = labels;
        initializeSlide();
    }

    public String getSelectedLabel() {
        return labels[selectedIndex];
    }

    public int getSelectedIndex() {
        return selectedIndex;
    }


    /************************************
     * Private methods
     ************************************/


    private void initializeSlide() {
        //The sliding pane
        slideWidth = viewWidth;
        slideHeight = UIC.gridY * labels.length;
        slideYMin = viewHeight - slideHeight;
        slide = context.createGraphics(viewWidth, slideHeight, PApplet.JAVA2D);
        slide.beginDraw();
        renderSlide();
    }


    private void renderList() {
        // cropped view onto the canvas
        PImage crop = slide.get(0, -slideY, viewWidth, viewHeight);
        context.image(crop, viewX, viewY);

        //shade variables
        int gradientA = context.color(0, 0, 0, 1);
        int gradientB = context.color(0, 0, 0, 64);

        //adding shade to top
        int topShadeY = -shadeHeight - slideY;
        if (topShadeY > 0) topShadeY = 0;
        for (int i = 0; i <= topShadeY + shadeHeight - 1; i++) {
            float inter = PApplet.map(i, topShadeY, topShadeY + shadeHeight, 0, 1);
            int c = context.lerpColor(gradientB, gradientA, inter);
            context.stroke(c);
            context.line(viewX, i + viewY, viewWidth + viewX, i + viewY);
        }

        //shading the bottom
        int minBottomY = viewHeight - shadeHeight;
        int bottomPortionHidden = slideHeight + slideY - viewHeight;
        int bottomShadeY = viewHeight - bottomPortionHidden;
        if (bottomShadeY < minBottomY) bottomShadeY = minBottomY;
        for (int i = bottomShadeY; i <= viewHeight - 1; i++) {
            float inter = PApplet.map(i, bottomShadeY, bottomShadeY + shadeHeight, 0, 1);
            int c = context.lerpColor(gradientA, gradientB, inter);
            context.stroke(c);
            context.line(viewX, i + viewY, viewWidth + viewX, i + viewY);
        }


        //drawing a top line that gets erased
        context.stroke(UIC.BG_LINE_COLOR);
        context.line(viewX, viewY, viewX + viewWidth, viewY);

        // shading right
        if (hasShadow) UIF.gradientX(viewX + viewWidth, viewY, UIC.gridXHalf, viewHeight, UIC.GRAD_B, UIC.GRAD_A);

    }


    private void renderSlide() {
        slide.background(UIC.BG_COLOR);
        slide.strokeWeight(1);
        slide.stroke(UIC.BG_LINE_COLOR);

        // vertical lines
        int xLoc = 0;
        while (xLoc < slideWidth) {
            slide.line(xLoc, 0, xLoc, slideHeight);
            xLoc += UIC.gridX;
        }

        // horizontal lines
        int yLoc = 0;
        while (yLoc < slideHeight) {
            slide.line(0, yLoc, slideWidth, yLoc);
            yLoc += UIC.gridY;
        }

        //rendering our labels
        slide.noStroke();
        slide.textFont(UIC.arial, UIC.fontSize);
        for (int index = 0; index < labels.length; index++) {
            if (index != selectedIndex) {
                //draw text
                String label = labels[index];
                slide.fill(UIC.TEXT_COLOR);
                UIF.drawTextGraphics(UIC.gridX, index * UIC.gridY, viewWidth, label, UIC.ALIGN_LEFT, slide);
            }
        }


        //selected button
        int highlightColor = UIF.getTranslucent(UIC.MAIN_COLOR, 0.8);
        slide.fill(highlightColor);
        slide.rect(UIC.gridX, UIC.gridY * selectedIndex, viewWidth - UIC.gridX, UIC.gridY);
        slide.fill(UIC.WHITE);
        UIF.drawTextGraphics(UIC.gridX, selectedIndex * UIC.gridY, viewWidth, labels[selectedIndex], UIC.ALIGN_LEFT, slide);
    }


    @SuppressWarnings("unused")
    private void drawButtonHilight() {
        context.fill(241, 94, 31, 180);
        context.rect(UIC.gridX * 4, (selectedIndex + 2) * UIC.gridY, UIC.gridX * 3, UIC.gridY);
        context.fill(241, 94, 31, 128);
        context.rect(UIC.gridX * 7, (selectedIndex + 2) * UIC.gridY, UIC.gridX * 1, UIC.gridY);

    }


    /************************************
     * Event handlers
     ************************************/


    public void mousePressed() {
        mousePressedY = context.mouseY;
        slideYOld = slideY;
    }


    public void mouseDragged() {
        if (labels.length * UIC.gridY > viewHeight) {
            slideY = slideYOld + (context.mouseY - mousePressedY);
            if (slideY > 0) slideY = 0;
            if (slideY < slideYMin) slideY = slideYMin;
        }
    }

    public void mouseReleased() {
        slideY = PApplet.round(slideY / UIC.gridY) * UIC.gridY;
        int newSelectedIndex = (int) (PApplet.floor((PApplet.abs(slideY) + (context.mouseY - viewY)) / UIC.gridY));
        if (newSelectedIndex >= labels.length) newSelectedIndex = labels.length - 1;
        if (PApplet.abs(context.mouseY - mousePressedY) < UIC.gridYFourth) {
            if (newSelectedIndex != selectedIndex) {
                selectedIndex = newSelectedIndex;
                listener.viewStateChanged(this);
                renderSlide();
            }

        }

    }

}
