/*
 * This represents one camphor boat. It contains a position, velocity,
 * and information about the surrounding boats.
 */
public class Boat {
	double R, Lrc;

	public double x;
	public double v;
	public double front, back;
	public int index, next, last;
	public boolean jammed;

	public Boat(double nR, double nLrc, int ind){
		R = nR;
		Lrc = nLrc;
		x = 0.0;
		v = 0.0;
		front = 0.0;
		back = 0.0;
		index = ind;
		next = 0;
		last = 0;
		jammed = false;
	}

	public Boat(int ind){
		R = 0.0;
		Lrc = 0.0;
		x = 0.0;
		v = 0.0;
		front = 0.0;
		back = 0.0;
		index = ind;
		next = 0;
		last = 0;
		jammed = false;
	}
	
	public void copyIn(Boat nb){
		x = nb.x;
		v = nb.v;
		
		front = nb.front;
		back = nb.back;
		
		index = nb.index;
		next = nb.next;
		last = nb.last;
	}
}
