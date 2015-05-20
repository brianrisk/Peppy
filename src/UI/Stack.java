package UI;

import processing.core.PApplet;

public class Stack extends Adapter implements Listener {
	
	int viewX = 0;
	int viewY = 0;
	PApplet context;
	
	ViewList list;
	View [] panes;
	
	public Stack(PApplet context) {
		this.context = context;
		list = new ViewList(context, this, Constants.gridX * 0, Constants.gridY * 1 + viewY, Constants.gridX * 4, Constants.gridY * 11);
		String [] menuLabels = {
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
		super.render();
		
		//put in the title
		context.fill(Constants.WHITE);

		//gradient colors
		int gradientA = context.color(0,0,0,1);
		int gradientB = context.color(0,0,0,51);
		
		//list gradient
		Functions.gradientX(list.viewX + list.viewWidth,list.viewY,Constants.gridXHalf,list.viewHeight,gradientB,gradientA);
		
		//menu gradient
		Functions.gradientY(0,Constants.gridYHalf + viewY,Constants.gridX * 16,Constants.gridYHalf,gradientA,gradientB);
		
	}
	

	
	
	private void renderBackground() {
		context.fill(Constants.BG_COLOR);
		context.rect(0,0,context.width, context.height);
			
		context.strokeWeight(1); 
		context.stroke(Constants.BG_LINE_COLOR);
	
		int xLoc = Constants.gridX;
		while (xLoc < context.width) {
			context.line(xLoc, 0, xLoc, context.height);
			xLoc += Constants.gridX;
		}
	
		int yLoc = Constants.gridY + viewY;
		while (yLoc < context.height) {
			context.line(0, yLoc, context.width, yLoc);
			yLoc += Constants.gridY;
		}
	}

}
