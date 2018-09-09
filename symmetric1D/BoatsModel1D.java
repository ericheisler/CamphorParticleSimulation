import javax.swing.*;
import java.util.*;

/**
 * This contains all of the model information and math
 * @author Heisler
 */
class BoatsModel1D {

	//These are parts of the model
	double R, m, vis, L, rc, D, volD, k, volk, a, vola, beta, csat, maxdgamma;
	int N;

	int nx;
	double dx;
	double[] xg;

	long nt, currentT;
	double T, dt;

	ArrayList<Boat1D> b;
	double[] c;
	double[] volc;
	double[] tmpcv;

	//for the modelEquations function
	double background;
	double[] supply;
	double[] Tf, Tb;
	int indb, indf;
	int rindb, rindf;
	double boffset, foffset;
	double rboffset, rfoffset;
	double avecb, avecf;

	// for the RK4 method
	double[][][] rk4Data;
	
	//for the crank nicolson method
	double upper, lower, diag;
	double bupper, blower, bdiag;
	double[] d;
	double kdt;
	double kvoldt;
	double adt;
	double[] fEven, fOdd;
	boolean evenStep;
	double[] mc, md, y, q;
	
	//for the volume crank nicolson method
	double vupper, vlower, vdiag;
	double vbupper, vblower, vbdiag;
	double[] vd;
	double vkdt;
	double vadt;
	double[] vfEven, vfOdd;
	boolean vevenStep;
	double[] vmc, vmd, vy, vq;
	
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
	public BoatsModel1D(){
		//here are some default parameters
		R = 20;
		m = 0.002;
		vis = 0.006;
		L = 0.3;
		rc = 0.15;
		D = 1;
		volD = 0.01;
		k = 0.1;
		volk = 0.1;
		a = 0;
		vola = 1;
		beta = 1;
		csat = 10;
		maxdgamma = 22;

		N = 0;
		T = 1;

		//discretization
		nx = (int)(50*R);
		dx = R/(nx-1);
		xg = new double[nx];
		for(int i=0; i<nx; i++){
			xg[i] = dx*i;
		}

		currentT = 0;
		evenStep = true;
		nt = (long)(T*50000);
		dt = T/(nt-1);

		// allocate all the necessary space
		b = new ArrayList<Boat1D>();
		c = new double[nx];
		volc = new double[nx];
		tmpcv = new double[nx];
		
		background = 0.0;
		Tf = new double[0];
		Tb = new double[0];

		rk4Data = new double[5][0][2];

		kdt = (k+volk)*dt/2;
		kvoldt = volk*dt/2;
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
		
		// volume version
		vkdt = volk*dt/2;
		vadt = vola*dt/2;
		vdiag = -volD*dt/dx/dx;
		vupper = volD*dt/2/dx/dx;
		vlower = volD*dt/2/dx/dx;
		vbupper = -vupper;
		vblower = -vlower;
		vbdiag = -vdiag+1;
		
		vd = new double[nx];
		vfEven = new double[nx];
		vfOdd = new double[nx];
		vmc = new double[nx];
		vmd = new double[nx];
		vy = new double[nx];
		vq = new double[nx];
		
		withNoise = false;
		noiseStrength = 0.0;
		noiseMaker = new Random();

		startind = 0;
		started = true;

		//initialize the arrays
		for(int i=0; i<nx; i++){
			c[i] = 0;
			volc[i] = 0;
			tmpcv[i] = 0.0;
		}

	}
	
	/**
     * This is the equation of motion
	 * y'=f(y)
	 * it writes the values of f(y)*dt directly to y
     * @param y
     */
    public double[][] motionEquation(double[][] xandv){
		
		// determine c at the front and back
		for (int i=0; i<N; i++){
			indb = (int)Math.floor((xandv[i][0]-L/2)/dx);		// find the indexes in xg for the back
			indf = (int)Math.floor((xandv[i][0]+L/2)/dx);		// and front of each boat
			boffset = (xandv[i][0]-L/2)/dx - indb;			//find the offsets
			foffset = (xandv[i][0]+L/2)/dx - indf;
			
			if(indb >= nx-1){ indb = indb - (nx-1); }	// for periodic boundary
			if(indf >= nx-1){ indf = indf - (nx-1); }
			if(indb < 0){ indb = indb + (nx-1); }		// for periodic boundary
			if(indf < 0){ indf = indf + (nx-1); }

			avecb = c[indb]*(1-boffset) + c[(indb+1)%nx]*boffset;
			avecf = c[indf]*(1-foffset) + c[(indf+1)%nx]*foffset;

			// calculate the surface tension
			Tf[i] = maxdgamma/(beta*beta*avecf*avecf + 1); //surface tension at front
			Tb[i] = maxdgamma/(beta*beta*avecb*avecb + 1); //surface tension at back
		}
		
		//finally calculate the diff. equations
		for (int i=0; i<N; i++){
			xandv[i][0] = xandv[i][1];
			xandv[i][1] = (-vis/m*xandv[i][1] + L/m*(Tf[i]-Tb[i]));	// for dv/dt
			// for dx/dt
		}

		return xandv;
	}
	
	// sets up the F vector
	public void setUpF(double[] f, double[][] pos){
		// clear the vector
		for (int i=0; i<f.length; i++){
			f[i] = 0.0;
		}
		
		//////////
		// no source on the surface
		// but this is used for the surface under the boat
		/////////
		
		// add constant sources
		for (int i=0; i<N; i++){
			rindb = (int)Math.floor((pos[i][0]-rc)/dx);			//find the camphor indices
			rindf = (int)Math.floor((pos[i][0]+rc)/dx);
			rboffset = (pos[i][0]-rc)/dx - rindb;				//find the offsets
			rfoffset = (pos[i][0]+rc)/dx - rindf;
			
			// add to the supply around each source
			int tmpj;
			for(int j=rindb+1; j<=rindf; j++){
				tmpj = j;
				if(j<0){ tmpj = j+(nx-1); }
				if(j>=nx-1){ tmpj = j-(nx-1); }
				f[tmpj] = 1;
			}
			//for the edges
			if(rindb<0){ rindb = rindb+(nx-1); }
			f[rindb] = (1-rboffset);			
			if(rindf>=nx-2){ rindf = rindf-(nx-1); }
			f[rindf+1] = rfoffset;
		}
		
	}
	
	// sets up the VF vector
	public void setUpVF(double[] f, double[][] pos){
		// clear the vector
		for (int i=0; i<f.length; i++){
			f[i] = 0.0;
		}
		// add constant sources
		for (int i=0; i<N; i++){
			rindb = (int)Math.floor((pos[i][0]-rc)/dx);			//find the camphor indices
			rindf = (int)Math.floor((pos[i][0]+rc)/dx);
			rboffset = (pos[i][0]-rc)/dx - rindb;				//find the offsets
			rfoffset = (pos[i][0]+rc)/dx - rindf;
			
			// add to the supply around each source
			int tmpj;
			for(int j=rindb+1; j<=rindf; j++){
				tmpj = j;
				if(j<0){ tmpj = j+(nx-1); }
				if(j>=nx-1){ tmpj = j-(nx-1); }
				f[tmpj] = (csat-volc[tmpj])*(csat-volc[tmpj])*(csat-volc[tmpj])/csat/csat/csat;
				//f[tmpj] = (csat-c[tmpj])*(csat-c[tmpj])*(csat-c[tmpj])/csat/csat/csat + k/a*c[tmpj];
			}
			//for the edges
			if(rindb<0){ rindb = rindb+(nx-1); }
			f[rindb] = (1-rboffset)*(csat-volc[rindb])*(csat-volc[rindb])*(csat-volc[rindb])/csat/csat/csat;			
			if(rindf>=nx-2){ rindf = rindf-(nx-1); }
			f[rindf+1] = rfoffset*(csat-volc[rindf+1])*(csat-volc[rindf+1])*(csat-volc[rindf+1])/csat/csat/csat;
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
	
	public void solveSystemVolume(){
		// uses sherman morrison method
		// solve Bc = d by By=d, Bq=u, c=y-(vy)/(1+vq)q
		// u = [-1..lower]
		// v = [1..-upper]    (maybe change some signs if it doesn't work)
		// NOTE: NOT using the end point nx-1. only consider 0..nx-2
		// for efficiency
		double id = 0.0;
		// solve for y
		vmc[0] = vbupper/(vbdiag+1);
		vmd[0] = vd[0]/(vbdiag+1);
		for (int i = 1; i < nx-2; i++){
			id = 1/(vbdiag - vmc[i-1]*vblower);
			vmc[i] = vbupper*id;
			vmd[i] = (vd[i] - vmd[i-1]*vblower)*id;
        }
		id = 1/(vbdiag + vbupper*vblower - vmc[nx-3]*vblower);
		vmc[nx-2] = vbupper*id;
		vmd[nx-2] = (vd[nx-2] - vmd[nx-3]*vblower)*id;
		
        vy[nx-2] = vmd[nx-2];
        for (int i=nx-3; i>=0; i--){
			vy[i]=vmd[i]-vmc[i]*vy[i+1];
		}
		
		//solve for q
		vmc[0] = vbupper/(vbdiag+1);
		vmd[0] = -1.0/(vbdiag+1);
		for (int i = 1; i < nx-2; i++){
			id = 1/(vbdiag - vmc[i-1]*vblower);
			vmc[i] = vbupper*id;
			vmd[i] = -vmd[i-1]*vblower*id;
        }
		id = 1/(vbdiag + vbupper*vblower - vmc[nx-3]*vblower);
		vmc[nx-2] = vbupper*id;
		vmd[nx-2] = (vblower - vmd[nx-3]*vblower)*id;
		
        vq[nx-2] = vmd[nx-2];
        for (int i=nx-3; i>=0; i--){
			vq[i]=vmd[i]-vmc[i]*vq[i+1];
		}
		
		//find the products
		double vvy = vy[0]-vy[nx-2]*vbupper;
		double vvq = vq[0]-vq[nx-2]*vbupper;
		
		//finally solve for c
		for (int i = 0; i < nx-1; i++){
			volc[i] = vy[i] - (vvy/(1+vvq))*vq[i];
        }
		volc[nx-1] = volc[0];
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// this is the HYBRID time integrator 
	// It does: 1. RK4 on motion eq.
	//			2. Crank-Nicolson on camphor eq.
	//			3. RK4 again with updated camphor
	// It is 4th order for motion, 2nd order for camphor
	// does that even make sense?
	// it does 1 time step per call
	//////////////////////////////////////////////////////////////////////////////
	public void timeStep(){
		// initial RK4 step
		rk4(false);
		
		// update F
		if(evenStep){ //step is even
			setUpF(fEven, rk4Data[0]);
			setUpVF(vfEven, rk4Data[0]);
		}else{
			setUpF(fOdd, rk4Data[0]);
			setUpVF(vfOdd, rk4Data[0]);
		}
		
		// crank-nicolson
		crankNicolson();
		
		// RK4 again with updated c
		rk4(true);
		
		// update F
		if(evenStep){ //step is even
			setUpF(fEven, rk4Data[0]);
			setUpVF(vfEven, rk4Data[0]);
		}else{
			setUpF(fOdd, rk4Data[0]);
			setUpVF(vfOdd, rk4Data[0]);
		}
		
		// increment time
		currentT++;
		evenStep = !evenStep;
	}
	
	// The RK4 integrator for the motion equations
	public void rk4(boolean complete){
		Boat1D bn;
		// perform a time step
		//prepare temp0 and step1
		for(int i=0; i<N; i++){
			bn = b.get(i);
			rk4Data[0][i][0] = bn.x;
			rk4Data[0][i][1] = bn.v;
			
			rk4Data[1][i][0] = bn.x;
			rk4Data[1][i][1] = bn.v;
		}
		rk4Data[1] = motionEquation(rk4Data[1]);	//first
		for(int i=0; i<N; i++){
			rk4Data[1][i][0] = rk4Data[1][i][0]*dt;
			rk4Data[1][i][1] = rk4Data[1][i][1]*dt;
		}
		
		for(int i=0; i<N; i++){
			rk4Data[2][i][0] = rk4Data[0][i][0] + rk4Data[1][i][0]/2;
			rk4Data[2][i][1] = rk4Data[0][i][1] + rk4Data[1][i][1]/2;
		}
		rk4Data[2] = motionEquation(rk4Data[2]);	//second
		for(int i=0; i<N; i++){
			rk4Data[2][i][0] = rk4Data[2][i][0]*dt;
			rk4Data[2][i][1] = rk4Data[2][i][1]*dt;
		}
		
		for(int i=0; i<N; i++){
			rk4Data[3][i][0] = rk4Data[0][i][0] + rk4Data[2][i][0]/2;
			rk4Data[3][i][1] = rk4Data[0][i][1] + rk4Data[2][i][1]/2;
		}
		rk4Data[3] = motionEquation(rk4Data[3]);	//third
		for(int i=0; i<N; i++){
			rk4Data[3][i][0] = rk4Data[3][i][0]*dt;
			rk4Data[3][i][1] = rk4Data[3][i][1]*dt;
		}
		
		for(int i=0; i<N; i++){
			rk4Data[4][i][0] = rk4Data[0][i][0] + rk4Data[3][i][0];
			rk4Data[4][i][1] = rk4Data[0][i][1] + rk4Data[3][i][1];
		}
		rk4Data[4] = motionEquation(rk4Data[4]);	//fourth
		for(int i=0; i<N; i++){
			rk4Data[4][i][0] = rk4Data[4][i][0]*dt;
			rk4Data[4][i][1] = rk4Data[4][i][1]*dt;
		}
		
		// store the step in temp0
		for(int i=0; i<N; i++){
			rk4Data[0][i][0] += (rk4Data[1][i][0] + 2*rk4Data[2][i][0] + 2*rk4Data[3][i][0] + rk4Data[4][i][0])/6;
			rk4Data[0][i][1] += (rk4Data[1][i][1] + 2*rk4Data[2][i][1] + 2*rk4Data[3][i][1] + rk4Data[4][i][1])/6;
		}
		
		// periodic boundary conditions
		for(int i=0; i<N; i++){
			if(rk4Data[0][i][0] >= R){
				rk4Data[0][i][0] -= R;
			}
			if(rk4Data[0][i][0] < 0){
				rk4Data[0][i][0] += R;
			}
		}
		
		
		if(complete){
			// if this is the final step to be stored
			for(int i=0; i<N; i++){
				bn = b.get(i);
				bn.x = rk4Data[0][i][0];
				bn.v = rk4Data[0][i][1];
			}
			
			// Take care of some annoying details
			for (int j=0; j<N; j++){
				update(j);
			}
			
			//first enforce the no-passing rule
			int bumpi = 0;
			for(Boat1D bnn: b){
				if(bnn.front > 0 && bnn.back > 0){
					bumpi = bnn.index; // this will be the index of the starting point
					break;
				}
			}
			for (int j=0; j<N; j++){
				if(b.get(bumpi).front < 0){
					bump(bumpi, -b.get(bumpi).front);
				}
				bumpi = (bumpi+1)%N;
			}
			
			//then enforce periodic boundary conditions again
			for(Boat1D bnn: b){
				if(bnn.x > R){
					bnn.x -= R;
					update(bnn.index);
					update(bnn.next);
					update(bnn.last);
				}
				if(bnn.x < 0){
					bnn.x += R;
					update(bnn.index);
					update(bnn.next);
					update(bnn.last);
				}
			}
		}
		
	}
	
	// The Crank-Nicolson integrator for the camphor equation
	public void crankNicolson(){
		////////////////////////////////////////////////////////
		// use the Crank-Nicolson method to step the camphor eq.
		////////////////////////////////////////////////////////
		//////////     Bc=d     ///////////
		// sruface:
		// d = [l*A-(k*dt/2-1)*I]*cn + a*dt/2*(Fnext + Fnow) + kv*dt/2*(tmpc+volc)
		// B = [-l*A+(k*dt/2+1)*I]
		// l = D*dt/2/dx^2
		// kdt = k*dt/2 - 1
		// adt = a*dt/2
		// k = ks+kv
		// 
		// bulk:
		// 
		////////////////////////////////////////////////////////
		// set up d NOTE: here A contains l
		
		// first solve the bulk equation
		// dc/dt = Dv d2c/dx2 + av Fv
		// store tmpcv for later
		for(int i=0; i<nx; i++){
			tmpcv[i] = volc[i];
		}
		vd[0] = vlower*volc[nx-2]+vdiag*volc[0]+vupper*volc[1] + volc[0] + vadt*(vfEven[0]+vfOdd[0]);
		vd[nx-1] = vd[0]; // WTF boundaries
		vd[nx-2] = vlower*volc[nx-3]+vdiag*volc[nx-2]+vupper*volc[0] + volc[nx-2] + vadt*(vfEven[nx-2]+vfOdd[nx-2]);
		for (int i=1; i<nx-2; i++){
			vd[i] = vlower*volc[i-1]+vdiag*volc[i]+vupper*volc[i+1] + volc[i] + vadt*(vfEven[i]+vfOdd[i]);
		}
		//B is already set, so solve the system for volc
		solveSystemVolume();
		
		// then solve the surface equation
		//
		d[0] = lower*c[nx-2]+diag*c[0]+upper*c[1] - (kdt-1)*c[0] + kvoldt*(tmpcv[0]+volc[0]);
		d[nx-1] = d[0]; // WTF boundaries
		d[nx-2] = lower*c[nx-3]+diag*c[nx-2]+upper*c[0] - (kdt-1)*c[nx-2] + kvoldt*(tmpcv[nx-2]+volc[nx-2]);
		for (int i=1; i<nx-2; i++){
			d[i] = lower*c[i-1]+diag*c[i]+upper*c[i+1] - (kdt-1)*c[i] + kvoldt*(tmpcv[i]+volc[i]);
		}
		//B is already set, so solve the system for c
		solveSystem();
		
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
			Boat1D newb = new Boat1D(0);
			newb.x = 0.0;
			newb.v = 0.0;
			newb.next = 0;
			newb.last = 0;
			b.add(0, newb);
			N++;
			placed = 0;
			update(0);
		}else {
			for(Boat1D bn: b){
				if(bn.front >= L){
					Boat1D newb = new Boat1D(bn.index+1);
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
				for(Boat1D bn: b){
					bn.index = i;
					bn.next = i+1;
					bn.last = i-1;
					i++;
				}
				b.get(0).last = N-1;
				b.get(N-1).next = 0;
				
				//enforce periodic boundaries
				for(Boat1D bn: b){
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
			Tf = new double[N];
			Tb = new double[N];
			double[][][] temprk4 = new double[5][N][2];
			for(int p=0; p<5; p++){
				for(int i=0; i<N-1; i++){
					for(int j=0; j<2; j++){
						temprk4[p][i][j] = rk4Data[p][i][j];
					}
				}
			}
			rk4Data = temprk4;
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
			for(Boat1D bn: b){
				bn.x = i*R/N;
				i++;
			}
			for(Boat1D bn: b){
				update(bn.index);
			}
		}
	}

	public void randomize(){
		if(N>0){
			Random rand = new Random();
			int i = 0;
			for(Boat1D bn: b){
				bn.x = i*R/N + rand.nextDouble()*(R/N-L-0.1);
				i++;
			}
			for(Boat1D bn: b){
				update(bn.index);
			}
		}
	}
	
	public int getN(){
		return N;
	}

	public ArrayList<Boat1D> getBoats(){
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
		nt = (long)(T*50000);
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
		volc = new double[nx];
		tmpcv = new double[nx];
		supply = new double[nx];
		double[][][] temprk4 = new double[5][N][2];

		kdt = (k+volk)*dt/2;
		kvoldt = volk*dt/2;
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
		
		// volume version
		vkdt = volk*dt/2;
		vadt = vola*dt/2;
		vdiag = -volD*dt/dx/dx;
		vupper = volD*dt/2/dx/dx;
		vlower = volD*dt/2/dx/dx;
		vbupper = -vupper;
		vblower = -vlower;
		vbdiag = -vdiag+1;
		
		vd = new double[nx];
		vfEven = new double[nx];
		vfOdd = new double[nx];
		vmc = new double[nx];
		vmd = new double[nx];
		vy = new double[nx];
		vq = new double[nx];
		
		//initialize c
		for(int i=0; i<nx; i++){
			c[i] = 0.0;
			volc[i] = 0.0;
			tmpcv[i] = 0.0;
		}
	}

	public double[] getParameters(){
		double[] params = new double[14];
		params[0] = R;
		params[1] = m;
		params[2] = vis;
		params[3] = L;
		params[4] = rc;
		params[5] = D;
		params[6] = volD;
		params[7] = k;
		params[8] = volk;
		params[9] = a;
		params[10] = vola;
		params[11] =  beta;
		params[12] = csat;
		params[13] = maxdgamma;
		return params;
	}

	public void setParameters(double[] params){
		if(R != params[0]){
			setR(params[0]);
		}
		m = params[1];
		vis = params[2];
		L = params[3];
		rc = params[4];
		D = params[5];
		volD = params[6];
		k = params[7];
		volk = params[8];
		a = params[9];
		vola = params[10];
		beta = params[11];
		csat = params[12];
		maxdgamma = params[13];
		// set some important stuff
		kdt = (k+volk)*dt/2;
		kvoldt = volk*dt/2;
		adt = a*dt/2;
		diag = -D*dt/dx/dx;
		upper = D*dt/2/dx/dx;
		lower = D*dt/2/dx/dx;
		bupper = -upper;
		blower = -lower;
		bdiag = -diag+kdt+1;
		// volume version
		vkdt = volk*dt/2;
		vadt = vola*dt/2;
		vdiag = -volD*dt/dx/dx;
		vupper = volD*dt/2/dx/dx;
		vlower = volD*dt/2/dx/dx;
		vbupper = -vupper;
		vblower = -vlower;
		vbdiag = -vdiag+1;
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
		evenStep = true;
		for(int i=N; i>0; i--){
			removeBoat();
		}
		for(int i=0; i<nx; i++){
			c[i] = 0.0;
			volc[i] = 0.0;
			tmpcv[i] = 0.0;
		}
	}
	
	public void changeBackground(double backc){
		background = backc;
	}
	
	public BoatsModel1D makeACopy(){
		BoatsModel1D newmodel = new BoatsModel1D();
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
	
	public void copyBoats(ArrayList<Boat1D> newb){
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
