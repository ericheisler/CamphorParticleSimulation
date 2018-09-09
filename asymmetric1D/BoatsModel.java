import javax.swing.*;
import java.util.*;

/**
 * This contains all of the model information and math
 * @author Heisler
 */
class BoatsModel {
	
	//These are parts of the model
	double R, m, vis, L, l, rc, D, k, a, beta;
	int N;
	
	int nx;
	double dx;
	double[] xg;
	
	long nt, currentT;
	double T, dt;
	
	ArrayList<Boat> b;
	double[] c;
	
	//for the modelEquations function
	double background;
	double[] cxx, supply;
	double[] Tf, Tb;
	int indb, indf;
	int rindb, rindf;
	double boffset, foffset;
	double rboffset, rfoffset;
	double avecb, avecf;
	
	// for the ABM4 method
	double[][] abm4tmp1, abm4tmp2, abm4tmp3, abm4tmp4, abm4tmp5;
	boolean abmPrepared;
	
	//for the crank nicolson method
	double upper, lower, diag;
	double bupper, blower, bdiag;
	double[] d;
	double kdt;
	double adt;
	double[] fEven, fOdd;
	boolean evenStep;
	double[] mc, md, y, q;
	
	// if we need some noise
	boolean withNoise;
	double noiseStrength;
	Random noiseMaker;
	
	int startind;
	boolean started;
	
	
	/**
	 * The constructor initializes all the model parts and prepares
	 * everything to run.
	 */
	public BoatsModel(){
		//here are some default parameters
		R = 45.5;
		m = 0.009;
		vis = 0.006;
		l = 0.6;
		L = 0.6;
		rc = 0.15;
		D = 1;
		k = 0.2;
		a = 20;
		beta = 0.2;
		
		N = 0;
		T = 1;
		
		//discretization
		nx = (int)(100*R);
		dx = R/(nx-1);
		xg = new double[nx];
		for(int i=0; i<nx; i++){
			xg[i] = dx*i;
		}
		
		currentT = 0;
		nt = (long)(T*500);
		dt = T/(nt-1);
		
		// allocate all the necessary space
		b = new ArrayList<Boat>();
		c = new double[nx];
		
		background = 0.0;
		Tf = new double[0];
		Tb = new double[0];
		
		abm4tmp1 = new double[2][0];
		abm4tmp2 = new double[2][0];
		abm4tmp3 = new double[2][0];
		abm4tmp4 = new double[2][0];
		abm4tmp5 = new double[2][0];
		abmPrepared = false;
		
		kdt = k*dt/2;
		adt = a*dt/2;
		diag = -D*dt/dx/dx;
		upper = D*dt/2/dx/dx;
		lower = D*dt/2/dx/dx;
		bupper = -upper;
		blower = -lower;
		bdiag = -diag+kdt+1;
		
		d = new double[nx];
		fEven = new double[nx];
		fOdd = new double[nx];
		mc = new double[nx];
		md = new double[nx];
		y = new double[nx];
		q = new double[nx];
		
		withNoise = false;
		noiseStrength = 0.0;
		noiseMaker = new Random();
		
		startind = 0;
		started = true;
		
		//initialize the arrays
		for(int i=0; i<nx; i++){
			c[i] = 0;
		}
	}
	
    /**
     * This is the equation of motion
	 * y'=f(y)
	 * it writes the values of f(y)*dt directly to y
     * @param y
     */
    public double[][] motionEquation(double[][] xandv){
		//// important ///////
		// v[i] = y[i]     //
		// x[i] = y[i+N]   //
		//////////////////////
		
		// determine c at the front and back
		for (int i=0; i<N; i++){
			indb = (int)Math.floor((xandv[0][i])/dx);		// find the indexes in xg for the back
			indf = (int)Math.floor((xandv[0][i]+L)/dx);		// and front of each boat
			boffset = (xandv[0][i])/dx - indb;			//find the offsets
			foffset = (xandv[0][i]+L)/dx - indf;
			
			if(indb >= nx-1){ indb = indb - (nx-1); }	// for periodic boundary
			if(indf >= nx-1){ indf = indf - (nx-1); }
			if(indb < 0){ indb = indb + (nx-1); }		// for periodic boundary
			if(indf < 0){ indf = indf + (nx-1); }
			
			avecb = c[indb]*(1-boffset) + c[(indb+1)%nx]*boffset;
			avecf = c[indf]*(1-foffset) + c[(indf+1)%nx]*foffset;
			
			// calculate the surface tension
			Tf[i] = 22/(beta*beta*avecf*avecf + 1); //surface tension at front
			Tb[i] = 22/(beta*beta*avecb*avecb + 1); //surface tension at back
		}
		
		//finally calculate the diff. equations
		for (int i=0; i<N; i++){
			xandv[0][i] = xandv[1][i];
			xandv[1][i] = (-vis/m*xandv[1][i] + l/m*(Tf[i]-Tb[i]));	// for dv/dt
			// for dx/dt
		}
		
		return xandv;
	}
	
	// sets up the F vector
	public void setUpF(double[] f, double[] pos){
		// clear the vector
		for (int i=0; i<f.length; i++){
			f[i] = 0.0;
		}
		// add constant sources
		for (int i=0; i<N; i++){
			rindb = (int)Math.floor((pos[i]-rc)/dx);			//find the camphor indices
			rindf = (int)Math.floor((pos[i]+rc)/dx);
			rboffset = (pos[i]-rc)/dx - rindb;				//find the offsets
			rfoffset = (pos[i]+rc)/dx - rindf;
			
			// add to the supply around each source
			int tmpj;
			for(int j=rindb+1; j<=rindf; j++){
				tmpj = j;
				if(j<0){ tmpj = j+(nx-1); }
				if(j>=nx-1){ tmpj = j-(nx-1); }
				f[tmpj] = 1.0;
			}
			//for the edges
			if(rindb<0){ rindb = rindb+(nx-1); }
			f[rindb] = (1-rboffset);			
			if(rindf>=nx-2){ rindf = rindf-(nx-1); }
			f[rindf+1] = rfoffset;
		}
	}
	
	public void solveSystem(){
		// uses sherman morrison method
		// solve Bc = d by By=d, Bq=u, c=y-(vy)/(1+vq)q
		// u = [-1..lower]
		// v = [1..-upper]    (maybe change some signs if it doesn't work)
		// NOTE: NOT using the end point nx-1. only consider 0..nx-2
		// for efficiency
		double id = 0.0;
		// solve for y
		mc[0] = bupper/(bdiag+1);
		md[0] = d[0]/(bdiag+1);
		for (int i = 1; i < nx-2; i++){
			id = 1/(bdiag - mc[i-1]*blower);
			mc[i] = bupper*id;
			md[i] = (d[i] - md[i-1]*blower)*id;
        }
		id = 1/(bdiag + bupper*blower - mc[nx-3]*blower);
		mc[nx-2] = bupper*id;
		md[nx-2] = (d[nx-2] - md[nx-3]*blower)*id;
		
        y[nx-2] = md[nx-2];
        for (int i=nx-3; i>=0; i--){
			y[i]=md[i]-mc[i]*y[i+1];
		}
		
		//solve for q
		mc[0] = bupper/(bdiag+1);
		md[0] = -1.0/(bdiag+1);
		for (int i = 1; i < nx-2; i++){
			id = 1/(bdiag - mc[i-1]*blower);
			mc[i] = bupper*id;
			md[i] = -md[i-1]*blower*id;
        }
		id = 1/(bdiag + bupper*blower - mc[nx-3]*blower);
		mc[nx-2] = bupper*id;
		md[nx-2] = (blower - md[nx-3]*blower)*id;
		
        q[nx-2] = md[nx-2];
        for (int i=nx-3; i>=0; i--){
			q[i]=md[i]-mc[i]*q[i+1];
		}
		
		//find the products
		double vy = y[0]-y[nx-2]*bupper;
		double vq = q[0]-q[nx-2]*bupper;
		
		//finally solve for c
		for (int i = 0; i < nx-1; i++){
			c[i] = y[i] - (vy/(1+vq))*q[i];
        }
		c[nx-1] = c[0];
	}
	
	// this is the HYBRID time integrator 
	// It does: 1. AB4 predictor on motion eq. only
	//			2. Crank-Nicolson on camphor eq.
	//			3. AM4 corrector on motion eq.
	// It is 4th order for motion, 2nd order for camphor
	// does that even make sense?
	// it does 1 time step per call
	public void timeStep(){
		/////////////// BE SURE TO COMPUTE THE FIRST 3 STEPS FIRST //////////////
		if(!abmPrepared){
			prepareABM4();
		}
		// now perform a time step
		
		////////////////////////////////////////////////////////
		// use AB4 on motion eq.
		////////////////////////////////////////////////////////
		//prepare temp1
		for (int i=0; i<N; i++){
			abm4tmp1[0][i] = b.get(i).x;
			abm4tmp1[1][i] = b.get(i).v;
		}
		
		abm4tmp1 = motionEquation(abm4tmp1);
		for (int i=0; i<N; i++){
			abm4tmp5[0][i] = b.get(i).x + (55*abm4tmp1[0][i] - 59*abm4tmp2[0][i] + 37*abm4tmp3[0][i] - 9*abm4tmp4[0][i])*dt/24;
			abm4tmp5[1][i] = b.get(i).v + (55*abm4tmp1[1][i] - 59*abm4tmp2[1][i] + 37*abm4tmp3[1][i] - 9*abm4tmp4[1][i])*dt/24 + noiseMaker.nextGaussian()*noiseStrength;
		}
		
		// now set up Feven or Fodd depending on the step number
		evenStep = currentT%2 == 0;
		if(evenStep){ //step is even
			setUpF(fEven, abm4tmp5[0]);
		}else{
			setUpF(fOdd, abm4tmp5[0]);
		}
		
		////////////////////////////////////////////////////////
		// use the Crank-Nicolson method to step the camphor eq.
		////////////////////////////////////////////////////////
		//////////     Bc=d     ///////////
		// d = [l*A-(k*dt/2-1)*I]*cn + a*dt/2*(Fnext + Fnow)
		// B = [-l*A+(k*dt/2+1)*I]
		// l = D*dt/2/dx^2
		// kdt = k*dt/2 - 1
		// adt = a*dt/2
		////////////////////////////////////////////////////////
		// set up d NOTE: here A contains l
		d[0] = lower*c[nx-2]+diag*c[0]+upper*c[1] - (kdt-1)*c[0] + adt*(fEven[0]+fOdd[0]);
		//d[nx-1] = lower*c[nx-2]+diag*c[nx-1]+upper*c[0] - (kdt-1)*c[nx-1] + adt*(fEven[nx-1]+fOdd[nx-1]);
		d[nx-1] = d[0]; // WTF boundaries
		d[nx-2] = lower*c[nx-3]+diag*c[nx-2]+upper*c[0] - (kdt-1)*c[nx-2] + adt*(fEven[nx-2]+fOdd[nx-2]);
		for (int i=1; i<nx-2; i++){
			d[i] = lower*c[i-1]+diag*c[i]+upper*c[i+1] - (kdt-1)*c[i] + adt*(fEven[i]+fOdd[i]);
		}
		//B is already set, so solve the system for c
		solveSystem();
		
		////////////////////////////////////////////////////////
		// use AM4 corrector on the motion eq.
		////////////////////////////////////////////////////////
		abm4tmp5 = motionEquation(abm4tmp5);
		// store values for the step
		for (int i=0; i<N; i++){
			b.get(i).x = b.get(i).x + (9*abm4tmp5[0][i] + 19*abm4tmp1[0][i] - 5*abm4tmp2[0][i] + abm4tmp3[0][i])*dt/24;
			b.get(i).v = b.get(i).v + (9*abm4tmp5[1][i] + 19*abm4tmp1[1][i] - 5*abm4tmp2[1][i] + abm4tmp3[1][i])*dt/24;
		}
		abm4tmp4 = abm4tmp3;
		abm4tmp3 = abm4tmp2;
		abm4tmp2 = abm4tmp1;
		
		// Take care of some annoying details
		for (int j=0; j<N; j++){
			update(j);
		}
		
		//first enforce the no-passing rule
		int bumpi = 0;
		for(Boat bn: b){
			if(bn.front > 0 && bn.back > 0){
				bumpi = bn.index; // this will be the index of the starting point
				break;
			}
		}
		for (int j=0; j<N; j++){
			if(b.get(bumpi).front < 0){
				bump(bumpi, -b.get(bumpi).front);
			}
			bumpi = (bumpi+1)%N;
		}
		
		//then enforce periodic boundary conditions
		for(Boat bn: b){
			if(bn.x > R){
				bn.x -= R;
				update(bn.index);
				update(bn.next);
				update(bn.last);
			}
			if(bn.x < 0){
				bn.x += R;
				update(bn.index);
				update(bn.next);
				update(bn.last);
			}
		}
		
		// and one more time, update F
		for (int i=0; i<N; i++){
			abm4tmp1[0][i] = b.get(i).x;
		}
		if(evenStep){ //step is even
			setUpF(fEven, abm4tmp1[0]);
		}else{
			setUpF(fOdd, abm4tmp1[0]);
		}
		
		//advance time
		currentT++;
	}
	
	// to prepare for ABM4
	// Just sets the first three steps to zero.
	// The boats are stationary and let's not advance time yet.
	public void prepareABM4(){
		//prepare temp4
		for (int i=0; i<N; i++){
			abm4tmp4[0][i] = b.get(i).x;
			abm4tmp4[1][i] = b.get(i).v;
		}
		//prepare temp3
		for (int i=0; i<N; i++){
			abm4tmp3[0][i] = b.get(i).x;
			abm4tmp3[1][i] = b.get(i).v;
		}
		//prepare temp2
		for (int i=0; i<N; i++){
			abm4tmp2[0][i] = b.get(i).x;
			abm4tmp2[1][i] = b.get(i).v;
		}
		abmPrepared = true;
	}
	
	//computes the distances between boats
	public void update(int index){
		int next = b.get(index).next;
		int last = b.get(index).last;
		if (N == 1){
			b.get(index).front = R-L;
			b.get(index).back = R-L;
		}else{
			if(b.get(next).x > b.get(index).x){
				b.get(index).front = b.get(next).x-b.get(index).x - L;
			}else{
				b.get(index).front = b.get(next).x-b.get(index).x + R - L;
			}if(b.get(last).x < b.get(index).x){
				b.get(index).back = b.get(index).x - b.get(last).x-L;
			}else{
				b.get(index).back = b.get(index).x - b.get(last).x-L+R;
			}
		}
	}
	
	//bumps the boats apart (recursively)
	public void bump(int index, double dist){
		int next = b.get(index).next;
		double frontspace = b.get(next).front;
		double backspace = b.get(index).back;
		double remaining;
		//if there is enough space, separate them and give average speed
		if(frontspace + backspace >= dist){
			if(frontspace < dist/2){
				remaining = Math.max(dist/2-frontspace, 0);
				b.get(next).x += Math.min(dist/2, frontspace);
				b.get(index).x -= remaining + dist/2;
			}else{
				remaining = Math.max(dist/2-backspace, 0);
				b.get(next).x += remaining + dist/2;
				b.get(index).x -= Math.min(dist/2, backspace);
			}
			b.get(index).v = (b.get(index).v+b.get(next).v)/2;
			b.get(next).v = b.get(index).v;
			update(index);
			update(next);
			update(b.get(next).next);
			update(b.get(index).last);
		}
		//if there is not space, bump forward and give average speed
		else{
			remaining = dist - (frontspace+backspace);
			
			//this gives slower speed
			//			if(b.get(index).v > b.get(next).v && b.get(next).v > 0){
			//				b.get(index).v = b.get(next).v;
			//			}
			//			if(b.get(index).v > b.get(next).v && b.get(next).v < 0){
			//				b.get(next).v = b.get(index).v;
			//			}
			
			//this gives avarage speed
			b.get(index).v = (b.get(index).v+b.get(next).v)/2;
			b.get(next).v = b.get(index).v;
			
			b.get(index).x -= backspace;
			b.get(index).front = 0.0;
			b.get(next).x += frontspace+remaining;
			b.get(next).back = 0.0;
			b.get(b.get(next).next).back = b.get(b.get(next).next).x - b.get(next).x;
			b.get(next).front = b.get(b.get(next).next).back;
			
			bump(next, remaining);
			
			//now reduce speed behind to prevent shockwaves
			//b.get(index).v = (b.get(index).v+b.get(next).v)/2;
			
			update(index);
			update(next);
			update(b.get(next).next);
			update(b.get(index).last);
		}
	}
	
	/**
	 * adds one boat to the system in the first available space
	 */
	public int addBoat(){
		int placed = -1;
		// stick it in the system
		if(N == 0){
			Boat newb = new Boat(0);
			newb.x = 0.0;
			newb.v = 0.0;
			newb.next = 0;
			newb.last = 0;
			b.add(0, newb);
			N++;
			placed = 0;
			update(0);
		}else {
			for(Boat bn: b){
				if(bn.front >= L){
					Boat newb = new Boat(bn.index+1);
					newb.x = bn.x+L;
					newb.v = bn.v;
					b.add(newb.index, newb);
					N++;
					placed = bn.index;
					break;
				}
			}
			if(placed >= 0){
				int i=0;
				for(Boat bn: b){
					bn.index = i;
					bn.next = i+1;
					bn.last = i-1;
					i++;
				}
				b.get(0).last = N-1;
				b.get(N-1).next = 0;
				
				//enforce periodic boundaries
				for(Boat bn: b){
					if(bn.x > R){
						bn.x -= R;
					}
					if(bn.x < 0){
						bn.x += R;
					}
				}
				
				update(placed);
				update(b.get(placed).next);
				update(b.get(placed).last);
			}
			if(placed < 0){
				JOptionPane.showMessageDialog(new JFrame(), "Couldn't add boat. Try again.");
			}
		}
		if(placed >= 0){
			// take care of some allocation details
			// keeping the old data.
			// this assumes that the added boat has been stationary.
			// I hope this doesn't fuck shit up
			Tf = new double[N];
			Tb = new double[N];
			abm4tmp1 = abm4tmp2;
			abm4tmp2 = new double[2][N];
			for(int i=0; i<N-1; i++){
				abm4tmp2[0][i] = abm4tmp1[0][i];
				abm4tmp2[1][i] = abm4tmp1[1][i];
			}
			abm4tmp1 = abm4tmp3;
			abm4tmp3 = new double[2][N];
			for(int i=0; i<N-1; i++){
				abm4tmp3[0][i] = abm4tmp1[0][i];
				abm4tmp3[1][i] = abm4tmp1[1][i];
			}
			abm4tmp1 = abm4tmp4;
			abm4tmp4 = new double[2][N];
			for(int i=0; i<N-1; i++){
				abm4tmp4[0][i] = abm4tmp1[0][i];
				abm4tmp4[1][i] = abm4tmp1[1][i];
			}
			abm4tmp1 = abm4tmp5;
			abm4tmp5 = new double[2][N];
			for(int i=0; i<N-1; i++){
				abm4tmp5[0][i] = abm4tmp1[0][i];
				abm4tmp5[1][i] = abm4tmp1[1][i];
			}
			abm4tmp1 = new double[2][N];
		}
		
		return N;
	}
	
	/**
	 * removes one boat from the system (from the last index)
	 */
	public int removeBoat(){
		if(N == 0){
			return N;
		}
		N--;
		b.remove(N);
		if(N > 1){
			b.get(N-1).next = 0;
			b.get(0).last = N-1;
			update(0);
			update(N-1);
		}
		
		// deallocate some stuff
		// no longer necessary
		
		return N;
	}
	
	public void homogenize(){
		if(N>0){
			int i = 0;
			for(Boat bn: b){
				bn.x = i*R/N;
				i++;
			}
			for(Boat bn: b){
				update(bn.index);
			}
		}
	}
	
	public void randomize(){
		if(N>0){
			Random rand = new Random();
			int i = 0;
			for(Boat bn: b){
				bn.x = i*R/N + rand.nextDouble()*(R/N-L-0.1);
				i++;
			}
			for(Boat bn: b){
				update(bn.index);
			}
		}
	}
	
	public int getN(){
		return N;
	}
	
	public ArrayList<Boat> getBoats(){
		return b;
	}
	
	public double[] getc(){
		return c;
	}
	
	public double getTime(){
		return currentT*dt;
	}
	
	public double getInterval(){
		return T;
	}
	
	public void setInterval(double newT){
		T = newT;
		nt = (long)(T*500);
	}
	
	public void setN(int newN){
		reset();
		for(int i=0; i<newN; i++){
			addBoat();
		}
	}
	
	public void setR(double newR){
		nx = (int)(100*newR);
		R = newR;
		dx = R/(nx-1);
		xg = new double[nx];
		for(int i=0; i<nx; i++){
			xg[i] = dx*i;
		}
		c = new double[nx];
		cxx = new double[nx];
		supply = new double[nx];
		abm4tmp1 = new double[2][N];
		abm4tmp2 = new double[2][N];
		abm4tmp3 = new double[2][N];
		abm4tmp4 = new double[2][N];
		abm4tmp5 = new double[2][N];
		abmPrepared = false;
		
		kdt = k*dt/2;
		adt = a*dt/2;
		diag = -D*dt/dx/dx;
		upper = D*dt/2/dx/dx;
		lower = D*dt/2/dx/dx;
		bupper = -upper;
		blower = -lower;
		bdiag = -diag+kdt+1;
		
		d = new double[nx];
		fEven = new double[nx];
		fOdd = new double[nx];
		mc = new double[nx];
		md = new double[nx];
		y = new double[nx];
		q = new double[nx];
		
		//initialize c
		for(int i=0; i<nx; i++){
			c[i] = 0.0;
		}
	}
	
	public double[] getParameters(){
		double[] params = new double[10];
		params[0] = R;
		params[1] = m;
		params[2] = vis;
		params[3] = l;
		params[4] = L;
		params[5] = rc;
		params[6] = D;
		params[7] = k;
		params[8] = a;
		params[9] =  beta;
		return params;
	}
	
	public void setParameters(double[] params){
		if(R != params[0]){
			setR(params[0]);
		}
		m = params[1];
		vis = params[2];
		l = params[3];
		L = params[4];
		rc = params[5];
		D = params[6];
		k = params[7];
		a = params[8];
		beta = params[9];
		// set some important stuff
		kdt = k*dt/2;
		adt = a*dt/2;
		diag = -D*dt/dx/dx;
		upper = D*dt/2/dx/dx;
		lower = D*dt/2/dx/dx;
		bupper = -upper;
		blower = -lower;
		bdiag = -diag+kdt+1;
	}
	
	public double[] getDisc(){
		double[] disc = new double[6];
		disc[0] = T;
		disc[1] = dt;
		disc[2] = (double)nt;
		disc[3] = R;
		disc[4] = dx;
		disc[5] = (double)nx;
		return disc;
	}
	
	public void reset(){
		currentT = 0;
		for(int i=N; i>0; i--){
			removeBoat();
		}
		for(int i=0; i<nx; i++){
			c[i] = 0.0;
		}
	}
	
	public void changeBackground(double backc){
		background = backc;
	}
	
	public BoatsModel makeACopy(){
		BoatsModel newmodel = new BoatsModel();
		newmodel.setParameters(getParameters());
		newmodel.setN(N);
		newmodel.setInterval(T);
		newmodel.copyc(c);
		newmodel.copyBoats(b);
		
		return newmodel;
	}
	
	public void copyc(double[] newc){
		if(c.length != newc.length){ 
			System.out.println("c was the wrong size?");
			return; 
		}
		System.arraycopy(newc,0,c,0,c.length);
	}
	
	public void copyBoats(ArrayList<Boat> newb){
		if(b.size() != newb.size()){ 
			System.out.println("b was the wrong size?");
			return; 
		}
		for(int i=0; i<b.size(); i++){
			b.get(i).copyIn(newb.get(i));
		}
	}
	
	public void moveABoat(int ind, double dist){
		// WARNING: may put a boat on top of another
		b.get(ind).x = (b.get(ind).x + dist)%R;
		
		update(ind);
		update(b.get(ind).next);
		update(b.get(ind).last);
	}
	
	public void kickABoat(int ind, double kick){
		b.get(ind).v = (b.get(ind).v + kick);
	}
	
	public void addNoise(boolean addit, double amount){
		withNoise = addit;
		noiseStrength = amount;
		noiseMaker = new Random();
	}
}
