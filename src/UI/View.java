package UI;

public abstract class View {
	
	public Listener listener;
	
	public int viewX;
	public int viewY;
	public int viewWidth;
	public int viewHeight;

	public boolean isActive = true;
	
	public abstract void render();
	public abstract void mousePressed();
	public abstract void mouseDragged();
	public abstract void mouseReleased();
	public abstract void updateLocation(int x, int y);
	
	/*
	 * returns if a point is inside the bounds of the component
	 */
	public boolean isIn(int pointX, int pointY) {
		return (pointX >= viewX && pointX < viewX + viewWidth) && (pointY >= viewY && pointY < viewY + viewHeight);
	}
	
	
}
