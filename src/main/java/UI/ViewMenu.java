package UI;

import processing.core.PApplet;

public class ViewMenu extends View {
	
	PApplet context;
	
	String [] items;
	
	int selectedIndex = 0;
	private int menuRight;
	private int menuWidth;
	
	/************************************
	 * Methods to override / implement
	 ************************************/
	 
	public void render() {
		//render background
		context.noStroke();
		context.fill(0,0,0,64);
		context.rect(viewX,viewY,menuWidth,viewHeight);
		
		//render items
		for (int index = 0; index < items.length; index++) {
			if (index == selectedIndex) {
				int highlightColor = UIF.getTranslucent(UIC.MAIN_COLOR, 0.8);
				context.fill(highlightColor);
				context.rect(viewX + (index * UIC.gridX * 2), viewY, UIC.gridX * 2, UIC.gridY);
				context.fill(UIC.WHITE);
			} else {
				context.fill(UIC.TEXT_COLOR);
			}
			UIF.drawText(viewX + (index * UIC.gridX * 2), viewY, UIC.gridX * 2, items[index], UIC.ALIGN_CENTER);
		}
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
	public ViewMenu(PApplet context, Listener listener, int x, int y, int width, int height) {
		this.context = context;
		this.listener = listener;
		this.viewX = x;
		this.viewY = y;
		this.viewWidth = width;
		this.viewHeight = height;
	}
	
	public void setItems(String [] items) {
		this.items = items;
		menuWidth = 2 * UIC.gridX * items.length;
		menuRight = viewX + menuWidth;
	}
	
	public String getItemLabel() {
		return items[selectedIndex];
	}
	
	public int getSelectedIndex() {
		return selectedIndex;
	}
	
	
	
	/************************************
	 * Event handlers
	 ************************************/
	 
	 
	public void mousePressed() {
		if (context.mouseX >= menuRight) return;
		selectedIndex = (int)((context.mouseX - viewX) / (UIC.gridX * 2));
		listener.viewStateChanged(this);
	}
	
	
	public void mouseDragged() {
		
	}
	
	public void mouseReleased() {

	}

}