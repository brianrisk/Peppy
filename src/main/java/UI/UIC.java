package UI;

import processing.core.PApplet;
import processing.core.PFont;
import Peppy.U;

/**
 * Constants
 * @author brianrisk
 *
 */
public class UIC {
	
	public static PApplet context;
	
	/* when a tablet */
	public static final int PLATFORM_COMPUTER = 0;
	public static final int PLATFORM_TABLET = 1;

	/* engage when testing */
	public static final boolean isTesting = false;

	/* DIMENSIONS */
	public static int gridX;
	public static int gridXHalf;
	public static int gridXFourth;
	public static int gridY;
	public static int gridYHalf;
	public static int gridYFourth;
	public static int gridMin;
	public static int gridMinHalf;
	public static int gridMinFourth;
	public static int fontSize;
	public static int fontSizeHalf;
	public static int fontSizeThird;

	/* FONTS */
	public static PFont arial;
	public static PFont arialBold;

	/* COLORS */
	public static int BG_COLOR;
	public static int BG_LINE_COLOR;
	public static int MAIN_COLOR;
	public static int WHITE;
	public static int BLACK;
	public static int TEXT_COLOR; //#6DA8CC
	public static int GRAD_A;
	public static int GRAD_B;

	/* ALIGNMENT */
	public static final int ALIGN_LEFT = 0;
	public static final int ALIGN_CENTER = 1;
	public static final int ALIGN_RIGHT = 2;

	//how many points in our pie circle
	public static final float wedgeRadianInc = PApplet.PI / 50;

	public static void init(PApplet context) {
		UIC.context = context;
		gridX = context.width / 16;
		gridXHalf = gridX / 2;
		gridXFourth = gridX / 4;
		gridY = context.height / 12;
		gridYHalf = gridY / 2;
		gridYFourth = gridY / 4;
		gridMin = PApplet.min(gridY, gridX);
		gridMinHalf = PApplet.min(gridYHalf, gridXHalf);
		gridMinFourth = PApplet.min(gridYFourth, gridXFourth);
		fontSize = gridYFourth;
		fontSizeHalf = fontSize / 2;
		fontSizeThird = fontSize / 3;
		
		/* FONTS */
		arial = context.createFont("Arial", fontSize);
		arialBold = context.createFont("Arial-BoldMT", fontSize);

		/* COLORS */
		BG_COLOR = context.color(31,39,47); //1f272F
		BG_LINE_COLOR = context.color(35, 51, 64);
		MAIN_COLOR = context.color(241,94,31);
		WHITE = context.color(255,255,255);
		BLACK = context.color (0,0,0);
		TEXT_COLOR = context.color(109,168,204); //#6DA8CC
		GRAD_A = context.color(0,0,0,1);
		GRAD_B = context.color(0,0,0,51);
	}

}
