package UI;

import processing.core.PApplet;

public class Stack extends Adapter implements Listener {

    int viewX = 0;
    int viewY = 0;
    PApplet context;

    ViewList list;
    View[] panes;

    public Stack(PApplet context) {
        this.context = context;
        list = new ViewList(context, this, UIC.gridX * 0, UIC.gridY * 1 + viewY, UIC.gridX * 4, UIC.gridY * 11);
        String[] menuLabels = {
                "Begin",
                "Sample",
                "Sequences",
                "Modifications",
                "Error",
                "Run"
        };
        list.setLabels(menuLabels);
        views.add(list);
    }

    public void viewStateChanged(View updatedView) {

    }

    public void render() {
        renderBackground();

        //put in the title
        context.fill(UIC.WHITE);
        context.textFont(UIC.arialBold);
        UIF.text(0, 0, 4, "PEPPY");
        context.textFont(UIC.arial);

        //menu gradient
        UIF.gradientY(0, UIC.gridYHalf + viewY, UIC.gridX * 16, UIC.gridYHalf, UIC.GRAD_A, UIC.GRAD_B);

        super.render();
    }


    private void renderBackground() {
        context.fill(UIC.BG_COLOR);
        context.rect(0, 0, context.width, context.height);

        context.strokeWeight(1);
        context.stroke(UIC.BG_LINE_COLOR);

        int xLoc = UIC.gridX;
        while (xLoc < context.width) {
            context.line(xLoc, 0, xLoc, context.height);
            xLoc += UIC.gridX;
        }

        int yLoc = UIC.gridY + viewY;
        while (yLoc < context.height) {
            context.line(0, yLoc, context.width, yLoc);
            yLoc += UIC.gridY;
        }
    }

}
