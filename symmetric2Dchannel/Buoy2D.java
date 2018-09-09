/*
 * This represents one camphor boat. It contains a position, velocity and index
 */
public class Buoy2D {

	public double x;
	public double y;
	public double vx;
	public double vy;
	
	public int index;

	public Buoy2D(int ind){
		x = 0.0;
		y = 0.0;
		vx = 0.0;
		vy = 0.0;
		
		index = ind;
	}
	
	// copies in the data from nb
	public void copyIn(Buoy2D nb){
		x = nb.x;
		y = nb.y;
		vx = nb.vx;
		vy = nb.vy;
		
		index = nb.index;
	}
}
