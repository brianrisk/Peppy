package UI;

import processing.core.PApplet;
import processing.core.PGraphics;

public class Functions {

	public static PApplet context;

	public static void init(PApplet context) {
		Functions.context = context;
	}

	/*
	 * Gradient
	 */
	public static void gradient(int x, int y, float w, float h, int c1, int c2, int axis ) {
		context.noFill();

		if (axis == 1) {  // Top to bottom gradient
			for (int i = y; i <= y+h - 1; i++) {
				float inter = PApplet.map(i, y, y+h, 0, 1);
				int c = context.lerpColor(c1, c2, inter);
				context.stroke(c);
				context.line(x, i, x+w, i);
			}
		}  
		else if (axis == 2) {  // Left to right gradient
			for (int i = x; i <= x+w - 1; i++) {
				float inter = PApplet.map(i, x, x+w, 0, 1);
				int c = context.lerpColor(c1, c2, inter);
				context.stroke(c);
				context.line(i, y, i, y+h);
			}
		}
	}

	public static void gradientX(int x, int y, float w, float h, int c1, int c2) {
		gradient(x,y,w,h,c1,c2,2);
	}

	public static void gradientY(int x, int y, float w, float h, int c1, int c2) {
		gradient(x,y,w,h,c1,c2,1);
	}

	/*
	 * PGraphics Gradient
	 */

	public static void gradientGraphics(int x, int y, float w, float h, int c1, int c2, int axis, PGraphics graphics ) {
		graphics.noFill();

		if (axis == 1) {  // Top to bottom gradient
			for (int i = y; i <= y+h - 1; i++) {
				float inter = PApplet.map(i, y, y+h, 0, 1);
				int c = context.lerpColor(c1, c2, inter);
				graphics.stroke(c);
				graphics.line(x, i, x+w, i);
			}
		}  
		else if (axis == 2) {  // Left to right gradient
			for (int i = x; i <= x+w - 1; i++) {
				float inter = PApplet.map(i, x, x+w, 0, 1);
				int c = context.lerpColor(c1, c2, inter);
				graphics.stroke(c);
				graphics.line(i, y, i, y+h);
			}
		}
	}

	public static void gradientXGraphics(int x, int y, float w, float h, int c1, int c2, PGraphics graphics) {
		gradientGraphics(x,y,w,h,c1,c2,2,graphics);
	}

	public static void gradientYGraphics(int x, int y, float w, float h, int c1, int c2, PGraphics graphics) {
		gradientGraphics(x,y,w,h,c1,c2,1,graphics);
	}


	/*
	 * Translucency
	 */

	public static int getTranslucent(int c, double percent) {
		if (percent > 1) percent = 1;
		return context.color(context.red(c), context.green(c), context.blue(c), (int) (percent * 255));
	}	


	/*
	 * Text
	 */

	public static void drawText(int xLoc, int yLoc, int pixelWidth, String label) {
		drawText(xLoc,xLoc,pixelWidth,label,Constants.ALIGN_LEFT);
	}

	public static void drawText(int xLoc, int yLoc, int pixelWidth, String label, int align) {
		
		int maxWidth = pixelWidth - Constants.gridXHalf;
		int labelWidth = (int)(context.textWidth(label));

		//printed is what we are actually printing, which may be truncated
		String printed = limitText(label, maxWidth);
		labelWidth = (int) (context.textWidth(printed));

		//draw text based on how it is aligned
		if (align == Constants.ALIGN_LEFT) {
			context.text(printed, xLoc + Constants.gridXFourth, yLoc + Constants.gridYHalf + Constants.fontSizeThird);
		}
		if (align == Constants.ALIGN_CENTER) {
			context.text(printed, xLoc + (pixelWidth - labelWidth) / 2, yLoc + Constants.gridYHalf + Constants.fontSizeThird);
		}
		if (align == Constants.ALIGN_RIGHT) {
			context.text(printed, xLoc + pixelWidth - (labelWidth + Constants.gridXFourth), yLoc + Constants.gridYHalf + Constants.fontSizeThird);
		}

	}

	public static void drawTextGraphics(int xLoc, int yLoc, int pixelWidth, String label, int align, PGraphics graphics) {
		int maxWidth = pixelWidth - Constants.gridXHalf;
		int labelWidth = (int)(context.textWidth(label));

		//printed is what we are actually printing, which may be truncated
		String printed = limitText(label, maxWidth);
		labelWidth = (int)(context.textWidth(printed));

		//draw text based on how it is aligned
		if (align == Constants.ALIGN_LEFT) {
			graphics.text(printed, xLoc + Constants.gridXFourth, yLoc + Constants.gridYHalf + Constants.fontSizeThird);
		}
		if (align == Constants.ALIGN_CENTER) {
			graphics.text(printed, xLoc + (pixelWidth - labelWidth) / 2, yLoc + Constants.gridYHalf + Constants.fontSizeThird);
		}
		if (align == Constants.ALIGN_RIGHT) {
			graphics.text(printed, xLoc + pixelWidth - (labelWidth + Constants.gridXFourth), yLoc + Constants.gridYHalf + Constants.fontSizeThird);
		}

	}

	public static String limitText(String label, int maxWidth) {
		String printed = label;
		int labelWidth = (int)(context.textWidth(label));

		//if the text is too wide, it is truncated at the maximum spot.
		if (labelWidth > maxWidth) {
			int elipsis = (int)(context.textWidth("..."));
			int desiredWidth = maxWidth - elipsis;
			int charIndex = label.length() - 2;
			while (labelWidth > desiredWidth) {
				printed = label.substring(0,charIndex--);
				labelWidth = (int)(context.textWidth(printed));
			}
			printed = printed + "...";
		}

		return printed;
	}


	/*
	 * Pie wedge for donut chart
	 */

	@SuppressWarnings("unused")
	private static void renderWedge(float middleX, float middleY, float radiusInner, float radiusOuter, float radiansA, float radiansB) {
		context.beginShape();
		float triangleWidth,triangleHeight;
		float radians;

		radians = radiansA;
		while (radians < radiansB) {
			triangleWidth = PApplet.cos(radians) *  radiusOuter;
			triangleHeight = PApplet.sin(radians) * radiusOuter;
			context.vertex(middleX + triangleWidth, middleY + triangleHeight);
			radians += Constants.wedgeRadianInc;
		}

		radians = radiansB;
		triangleWidth = PApplet.cos(radians) *  radiusOuter;
		triangleHeight = PApplet.sin(radians) * radiusOuter;
		context.vertex(middleX + triangleWidth, middleY + triangleHeight);

		while (radians > radiansA) {
			triangleWidth = PApplet.cos(radians) *  radiusInner;
			triangleHeight = PApplet.sin(radians) * radiusInner;
			context.vertex(middleX + triangleWidth, middleY + triangleHeight);
			radians -= Constants.wedgeRadianInc;
		}

		radians = radiansA;
		triangleWidth = PApplet.cos(radians) *  radiusInner;
		triangleHeight = PApplet.sin(radians) * radiusInner;
		context.vertex(middleX + triangleWidth, middleY + triangleHeight);

		context.endShape(PApplet.CLOSE);
	}


}
