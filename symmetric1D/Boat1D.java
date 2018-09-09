/*
 * This represents one camphor boat. It contains a position, velocity,
 * and information about the surrounding boats.
 */
public class Boat1D {

	public double x;
	public double v;
	public double front, back;
	public int index, next, last;
	public boolean jammed;

	public Boat1D(int ind){
		x = 0.0;
		v = 0.0;
		
		front = 0.0;
		back = 0.0;
		
		index = ind;
		next = 0;
		last = 0;
	}

	public void copyIn(Boat1D nb){
		x = nb.x;
		v = nb.v;
		
		front = nb.front;
		back = nb.back;
		
		index = nb.index;
		next = nb.next;
		last = nb.last;
	}
}
