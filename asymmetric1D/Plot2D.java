import java.awt.*;
import java.awt.geom.*;
import java.awt.image.*;
import java.util.*;

/**
 * An Image(BufferedImage) which represents a 2-D data plot
 *
 * @author H
 */
public class Plot2D extends BufferedImage {
	public int width;
	public int height;
	double[] xdata;
	double[][] ydata;
	double[][] loosePoints;
	double xmin, xmax, ymin, ymax;
	int[] pointType;
	int looseType;
	Vector<Color> colors;
	Color looseColor;
	String title, xlabel, ylabel;
	int topPad, xPad, yPad;

	boolean boxed, showAxes, xlimSet, ylimSet, lines;

	public static int DOT = 0;
	public static int CIRCLE = 1;
	public static int RING = 2;
	public static int SQUARE = 3;
	public static int BOX = 4;

	public Plot2D(int wid, int hig){
		super(wid, hig, BufferedImage.TYPE_INT_ARGB);
		width = wid;
		height = hig;
		colors = new Vector<Color>();
		boxed = true;
		showAxes = false;
		xlimSet = false;
		ylimSet = false;
		lines = false;
		topPad = 0;
		xPad = 0;
		yPad = 0;
	}

	public boolean addX(double[] newx){
		if(newx == null){
			return false;
		}
		xdata = newx;
		return true;
	}

	public boolean addY(double[] newy){
		if(newy == null){
			return false;
		}
		if(ydata == null){
			ydata = new double[1][newy.length];
			for(int j=0; j<newy.length; j++){
				ydata[0][j] = newy[j];
			}
			pointType = new int[1];
			pointType[0] = CIRCLE;
			colors.add(Color.RED);
			return true;
		}
		double[][] tmp = new double[ydata.length+1][ydata[0].length];
		for(int i=0; i<ydata.length; i++){
			for(int j=0; j<ydata[0].length; j++){
				tmp[i][j] = ydata[i][j];
			}
		}
		for(int j=0; j<newy.length; j++){
			tmp[tmp.length-1][j] = newy[j];
		}
		ydata = tmp;
		int[] tmp2 = new int[pointType.length+1];
		for(int i=0; i<pointType.length; i++){
			tmp2[i] = pointType[i];
		}
		tmp2[pointType.length] = CIRCLE;
		pointType = tmp2;
		colors.add(Color.RED);
		return true;
	}

	public void addY(double[] newy, int type, Color col){
		if(addY(newy)){
			pointType[pointType.length-1] = type;
			colors.set(colors.size()-1, col);
		}
	}
	
	public boolean addPoint(double px, double py){
		if(loosePoints == null){
			loosePoints = new double[1][2];
			loosePoints[0][0] = px;
			loosePoints[0][1] = py;
			looseType = CIRCLE;
			looseColor = Color.BLUE;
			return true;
		}
		int newind = loosePoints.length;
		double[][] tmp = loosePoints;
		loosePoints = new double[newind+1][2];
		loosePoints[newind][0] = px;
		loosePoints[newind][1] = py;
		System.arraycopy(tmp, 0, loosePoints, 0, tmp.length);
		tmp = null;
		return true;
	}
	
	public void addPoint(double px, double py, int type, Color col){
		if(addPoint(px, py)){
			looseType = type;
			looseColor = col;
		}
	}

	public void setXLimits(double min, double max){
		xmin = min;
		xmax = max;
		xlimSet = true;
	}

	public void setYLimits(double min, double max){
		ymin = min;
		ymax = max;
		ylimSet = true;
	}
	
	public void setTitle(String s){
		title = s;
		topPad = 15;
	}
	
	public void setLabels(String xs, String ys){
		if(xs != null){
			xlabel = xs;
			xPad = 20;
		}
		if(ys != null){
			ylabel = ys;
			yPad = 20;
		}
	}
	
	public void showAxes(boolean b){
		showAxes = b;
	}
	
	public void boxed(boolean b){
		boxed = b;
	}

	public void lines(boolean b){
		lines = b;
	}

	public void clearData(){
		xdata = null;
		ydata = null;
		loosePoints = null;
		pointType = null;
		colors.clear();
	}

	public double[] getExtremeValues() {
		if(xdata == null || ydata == null){
			return new double[] {0, 1, 0, 1};
		}
		double[] values = new double[4];
		if(!xlimSet){
			xmin = Double.MAX_VALUE-1;
			xmax = -xmin;
			for(int i = 0; i < xdata.length; i++) {
				xmax = Math.max(xmax, xdata[i]);
				xmin = Math.min(xmin, xdata[i]);
			}
			if(xmin >= xmax){
				xmin = Math.min(0.0, xmax-0.1);
			}
			values[0] = xmin;
			values[1] = xmax;
			//add a little something
			if(xmax - xmin > 1){
				xmax = Math.ceil(xmax);
				xmin = Math.floor(xmin);
			}else if(xmax - xmin > 0.1){
				xmax = Math.ceil(10*xmax)/10;
				xmin = Math.floor(10*xmin)/10;
			}
		}
        if(!ylimSet){
			ymin = Double.MAX_VALUE-1;
			ymax = -ymax;
			for(int i = 0; i < ydata.length; i++) {
				for(int j=0; j<ydata[i].length; j++){
					ymax = Math.max(ymax, ydata[i][j]);
					ymin = Math.min(ymin, ydata[i][j]);
				}
			}
			if(ymin >= ymax){
				ymin = Math.min(0.0, ymax-0.1);
			}
			values[2] = ymin;
			values[3] = ymax;
			//add a little something
			if(ymax - ymin > 1){
				ymax = Math.ceil(ymax);
				ymin = Math.floor(ymin);
			}else if(ymax-ymin > 0.1){
				ymax = Math.ceil(2*ymax)/2;
				ymin = Math.floor(2*ymin)/2;
			}
		}
		return values;
    }

	public void paint(){
		Graphics2D g = createGraphics();
		g.setColor(Color.WHITE);
		g.fillRect(0, 0, width, height);
		g.setColor(Color.BLACK);
		if(boxed){
			g.drawRect(0, 0, width-1, height-1);
		}
		if(xdata == null || ydata == null){
			g.drawString("empty plot", 10, (int)(height/2));
			return;
		}
		if(ydata.length == 0 || xdata.length == 0){
			g.drawString("empty plot", 10, (int)(height/2));
			return;
		}

		getExtremeValues();
		int pwidth = width-yPad;
		int pheight = height-topPad-xPad;
		if(showAxes){
			int xline = 0+xPad; //the x axis
			int yline = 0+yPad; //the y axis
			if(xmin < 0){
				yline = (int)((0.0-xmin)/(xmax-xmin)*pwidth+yPad);
				g.drawLine(yline, 0+topPad, yline, height-xPad); // the zero line y
			}
			if(ymin < 0){
				xline = (int)(topPad+pheight - (0.0-ymin)/(ymax-ymin)*pheight);
				g.drawLine(0+yPad, xline, width, xline); // the zero line x
			}
			g.drawLine(yPad, 0+topPad, yPad, height-xPad); //the box line y
			g.drawLine(0+yPad, height-xPad, width, height-xPad); //the box line x
			//draw tick marks
			g.drawLine(yline-5, (int)(pheight/4)+topPad, yline+5, (int)(pheight/4)+topPad);
			g.drawLine(yline-5, (int)(pheight/2)+topPad, yline+5, (int)(pheight/2)+topPad);
			g.drawLine(yline-5, (int)(pheight*3/4)+topPad, yline+5, (int)(pheight*3/4)+topPad);
			g.drawLine((int)(pwidth/4)+yPad, height-xPad-5, (int)(pwidth/4)+yPad, height-xPad+5);
			g.drawLine((int)(pwidth/2)+yPad, height-xPad-5, (int)(pwidth/2)+yPad, height-xPad+5);
			g.drawLine((int)(pwidth*3/4)+yPad, height-xPad-5, (int)(pwidth*3/4)+yPad, height-xPad+5);
			//draw numbers
			// the x axis numbers
			if(xlabel != null){
				String s = String.format("%.1f", xmin);
				g.drawString(s, yPad, height-2);
				s = String.format("%.1f", xmax);
				g.drawString(s, width - 10*s.length(), height-2);
			}
			
			//the y axis numbers
			if(ylabel != null){
				Font font = g.getFont();
				AffineTransform fat = new AffineTransform();
				fat.rotate(-Math.PI/2);
				Font rfont = font.deriveFont(fat);
				g.setFont(rfont);

				String s = String.format("%.3f", ymax);
				g.drawString(s, 15, 10*s.length());
				s = String.format("%.3f", ymin);
				g.drawString(s, 15, height - xPad);
				g.setFont(font);
			}
		}
		//draw labels
		if(title != null){
			int ts = (int)Math.max(0, (width-5*title.length())/2);
			g.drawString(title, ts, 12);
		}
		if(xlabel != null){
			int ts = (int)Math.max(0, (width-5*xlabel.length())/2);
			g.drawString(xlabel, ts, height-2);
		}
		if(ylabel != null){
			Font font = g.getFont();
			AffineTransform fat = new AffineTransform();
			fat.rotate(-Math.PI/2);
			Font rfont = font.deriveFont(fat);
			g.setFont(rfont);
			int ts = (int)Math.max(0, (height+5*ylabel.length())/2);
			g.drawString(ylabel, 12, ts);
			g.setFont(font);
		}


		//finally, draw the data
		int xpix, lxp;
		int[] ypix = new int[ydata.length];
		int[] lyp = new int[ydata.length];
		double lxd = 0.0;
		lxp = 0;
		lyp[0] = 0;
		for(int i=0; i<xdata.length; i++){
			xpix = (int)((xdata[i]-xmin)/(xmax-xmin)*pwidth+yPad);
			for(int j=0; j<ydata.length; j++){
				ypix[j] = (int)(topPad+pheight - (ydata[j][i]-ymin)/(ymax-ymin)*pheight);
				g.setColor(colors.get(j));
				if(pointType[j] == DOT){
					g.fillOval(xpix-1, ypix[j]-1, 3, 3);
				}else if(pointType[j] == CIRCLE){
					g.fillOval(xpix-2, ypix[j]-2, 5, 5);
				}else if(pointType[j] == RING){
					g.drawOval(xpix-3, ypix[j]-3, 7, 7);
				}else if(pointType[j] == SQUARE){
					g.fillRect(xpix-2, ypix[j]-2, 5, 5);
				}else if(pointType[j] == BOX){
					g.drawRect(xpix-3, ypix[j]-3, 7, 7);
				}
				if(lines && i>0 && xdata[i] > lxd){
					g.drawLine(lxp, lyp[j], xpix, ypix[j]);
				}
				lyp[j] = ypix[j];
			}
			lxd = xdata[i];
			lxp = xpix;
		}
		if(lines && xdata[0] > xdata[xdata.length-1]){
			xpix = (int)((xdata[0]-xmin)/(xmax-xmin)*pwidth+yPad);
			for(int j=0; j<ydata.length; j++){
				ypix[j] = (int)(topPad+pheight - (ydata[j][0]-ymin)/(ymax-ymin)*pheight);
				g.drawLine(lxp, lyp[j], xpix, ypix[j]);
			}
		}
		//and draw the loose points scaling them to the displayed range
		if(loosePoints != null){
			double lpmin = Double.MAX_VALUE-1;
			double lpmax = -lpmin;
			for(int i=0; i<loosePoints.length; i++){
				lpmax = Math.max(lpmax, loosePoints[i][1]);
				lpmin = Math.min(lpmin, loosePoints[i][1]);
			}
			for(int i=0; i<loosePoints.length; i++){
				xpix = (int)((loosePoints[i][0]-xmin)/(xmax-xmin)*pwidth+yPad);
				ypix[0] = (int)(topPad+pheight - (loosePoints[i][1]-lpmin)/(lpmax-lpmin)*pheight);
				g.setColor(looseColor);
				if(looseType == DOT){
					g.fillOval(xpix-1, ypix[0]-1, 3, 3);
				}else if(looseType == CIRCLE){
					g.fillOval(xpix-2, ypix[0]-2, 5, 5);
				}else if(looseType == RING){
					g.drawOval(xpix-3, ypix[0]-3, 7, 7);
				}else if(looseType == SQUARE){
					g.fillRect(xpix-2, ypix[0]-2, 5, 5);
				}else if(looseType == BOX){
					g.drawRect(xpix-3, ypix[0]-3, 7, 7);
				}
			}
		}
	}
}
