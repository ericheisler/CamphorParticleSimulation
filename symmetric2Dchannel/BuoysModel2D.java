import javax.swing.*;
import java.util.*;

/**
 * This contains all of the model information and math
 * @author Heisler
 */
class BuoysModel2D {

	//These are parts of the model
	public double R, W, m, vis, L, rc, D, k, a, beta, csat, maxdgamma;
	public double restitution;
	public int N;

	public int nx;			// x is long, y is short
	public int ny;
	public double dx; // dx and dy are same
	public double[] xg;
	public double[] yg;
	public int stepsperrc, stepspersec;

	public long nt, currentT;
	public double T, dt;

	public ArrayList<Buoy2D> b;
	public double[][] c;
	public double[][] distancesq;

	//for the modelEquations function
	double background;
	double[] boatEdgeST; // tension at 32 points evenly spaced around the edge of the boat
	double edgeAngle;
	double[] cosAngle, sinAngle;
	double edgex, edgey, deltax, deltay;
	public double[] netx, nety;
	
	// for the RK4 method
	double[][][] rk4Data;
	
	// for the ADI method
	double upper, lower, diag;
	double tmpm, tmpdiag;
	double lambda;
	double usefulnumber, adt;
	double [][]tmpc;
	double [] rhs;
	double [][] fn, fnh, fnp;
	
	double[] mc, md, y, q;
	
	// if we need some noise
	boolean withNoise;
	double noiseStrength;
	Random noiseMaker;
	
	// for statistics
	public boolean computeStats;
	public double[] initialx, initialy;
	public double[] meanDisp, meansqDisp;
	public double[] aveminDist;
	public double currentMinDist;
	public double[] avev;

	int startind;
	boolean started;


	/**
	 * The constructor initializes all the model parts and prepares
	 * everything to run.
	 */
	public BuoysModel2D(){
		//here are some default parameters
		R = 20;
		W = 0.4;
		m = 0.003;
		vis = 0.1;
		L = 0.3;
		rc = 0.15;
		D = 0.1;
		k = 0.1;
		a = 10;
		beta = 0.3333333333333;
		restitution = 0.5;
		csat = 10;
		maxdgamma = 22;

		N = 0;
		T = 1;

		//discretization
		stepsperrc = 5;
		stepspersec = 3000;
		nx = (int)(stepsperrc*R/rc); // gives 10 points across the camphor
		ny = (int)(stepsperrc*W/rc);
		dx = R/(nx-1);
		xg = new double[nx];
		yg = new double[ny];
		for(int i=0; i<nx; i++){
			xg[i] = dx*i;
		}
		for(int i=0; i<ny; i++){
			yg[i] = dx*i;
		}

		currentT = 0;
		nt = (long)(T*stepspersec);
		dt = T/(nt-1);

		// allocate all the necessary space
		b = new ArrayList<Buoy2D>();
		c = new double[nx][ny];
		
		background = 0.0;
		boatEdgeST = new double[64];
		edgeAngle = 3.14159265358979323846*2/64;
		cosAngle = new double[64];
		sinAngle = new double[64];
		for(int i=0; i<64; i++){
			cosAngle[i] = Math.cos(i*edgeAngle);
			sinAngle[i] = Math.sin(i*edgeAngle);
		}
		netx = new double[0];
		nety = new double[0];
		
		rk4Data = new double[5][4][0];
		
		distancesq = new double[0][0];
		
		tmpc = new double[nx][ny];
		rhs = new double[nx];
		fn = new double[nx][ny];
		fnh = new double[nx][ny];
		fnp = new double[nx][ny];
		
		mc = new double[nx];
		md = new double[nx];
		y = new double[nx];
		q = new double[nx];

		lambda = D*dt/2/dx/dx;
		usefulnumber = 1-2*lambda - dt*k/4;
		adt = dt/4;
		diag = 1.0+2*lambda + dt*k/4;
		upper = -lambda;
		lower = -lambda;
		tmpm = lower/diag;
		tmpdiag = diag - tmpm*upper;
		
		withNoise = true;
		noiseStrength = 0.1;
		noiseMaker = new Random();
		
		computeStats = false;
		initialx = new double[0];
		initialy = new double[0];
		meanDisp = new double[0];
		meansqDisp = new double[0];
		aveminDist = new double[0];
		currentMinDist = 0;
		avev = new double[0];

		startind = 0;
		started = true;

		//zero the arrays
		for(int i=0; i<nx; i++){
			for(int j=0; j<ny; j++){
				c[i][j] = 0;
			}
		}

	}
	
	/**
     * This is the equation of motion
	 * y'=f(y)
	 * it writes the values of f(y)*dt directly to y
     * @param y
     */
    public double[][] motionEquation(double[][] xyv){
		int blx, bly;
		
		// compute the net surface tension on each boat
		// then compute the RHS of the differential eq.
		for(Buoy2D bn: b){
			// use bilinear interpolation at 32 evenly spaced points 
			// around the edge of the boat
			netx[bn.index] = 0.0;
			nety[bn.index] = 0.0;
			for(int i=0; i<64; i++){
				edgex = bn.x + L/2*cosAngle[i];
				edgey = bn.y + L/2*sinAngle[i];
				if(edgex < 0){ edgex = 0; }
				if(edgey < 0){ edgey = 0; }
				if(edgex > R){ edgex = R; }
				if(edgey > W){ edgey = W; }
				
				blx = (int)(edgex/dx);
				bly = (int)(edgey/dx);
				if(blx > nx-2){ blx = nx-2; }
				if(bly > ny-2){ bly = ny-2; }
				
				deltax = edgex - blx*dx;
				deltay = edgey - bly*dx;
				
				boatEdgeST[i] = ( c[blx][bly]*(dx-deltax)*(dx-deltay) + c[blx+1][bly]*deltax*(dx-deltay) + c[blx][bly+1]*(dx-deltax)*deltay + c[blx+1][bly+1]*deltax*deltay )/dx/dx;
				
				boatEdgeST[i] = maxdgamma/(beta*beta*boatEdgeST[i]*boatEdgeST[i] + 1);
				
				netx[bn.index] += boatEdgeST[i]*L*edgeAngle/2*cosAngle[i];
				nety[bn.index] += boatEdgeST[i]*L*edgeAngle/2*sinAngle[i];
			}
			
			xyv[bn.index][0] = xyv[bn.index][2];
			xyv[bn.index][1] = xyv[bn.index][3];
			xyv[bn.index][2] = -vis/m*xyv[bn.index][2] + netx[bn.index]/m;
			xyv[bn.index][3] = -vis/m*xyv[bn.index][3] + nety[bn.index]/m;
		}
		
		return xyv;
	}
	
	// computes the distances between boats 
	public void computeDistances(){
		Buoy2D b1, b2;
		for(int i=0; i<N; i++){
			b1 = b.get(i);
			for(int j=i+1; j<N; j++){
				b2 = b.get(j);
				distancesq[i][j] = (b1.x-b2.x)*(b1.x-b2.x) + (b1.y-b2.y)*(b1.y-b2.y);
				distancesq[j][i] = distancesq[i][j];
			}
		}
	}
	
	// checks for and handles collisions
	public void doCollisions(){
		Buoy2D b1, b2;
		boolean somethingCollided = true;
		int tries = 0;
		
		while(somethingCollided && tries < 1){
			somethingCollided = false;
			tries++;
			
			for(int i=0; i<N; i++){
				b1 = b.get(i);
				// check for wall collisions
				// manually move the boats out of the wall
				if(b1.y < L/2){
					b1.y = L/2 + (L/2-b1.y)*restitution;
					if(b1.vy < 0){
						b1.vy = -b1.vy*restitution;
					}
					somethingCollided = true;
				}else if(b1.y > W-L/2){
					b1.y = W-L/2 - (b1.y-(W-L/2))*restitution;
					if(b1.vy > 0){
						b1.vy = -b1.vy*restitution;
					}
					somethingCollided = true;
				}
				if(b1.x < 0){
					b1.x += R;
					somethingCollided = true;
				}else if(b1.x >= R){
					b1.x -= R;
					somethingCollided = true;
				}
				
				// check for collisions between boats
				for(int j=i+1; j<N; j++){
					b2 = b.get(j);
					
					if(distancesq[i][j] < L*L){
						boatCollision(b1, b2);
						
						somethingCollided = true;
					}
				}
			}
			
			// recompute distances if needed
			if(somethingCollided && computeStats){ computeDistances(); }
		}
		
		if(somethingCollided){
			//computeDistances();
			//System.out.println("Warning: Couldn't handle collisions.");
		}
	}
	
	// handles a collision between boat1 and boat2
	// it is elastic
	public void boatCollision(Buoy2D b1, Buoy2D b2){
		// distance
		double overlap = L - Math.sqrt(distancesq[b1.index][b2.index]);
		/*
		// velocity in the direction of (dx, dy)
		double cv1 = (b1.vx*(b2.x-b1.x) + b1.vy*(b2.y-b1.y))/dist;
		double cv2 = (b2.vx*(b2.x-b1.x) + b2.vy*(b2.y-b1.y))/dist;
		
		// how long since collision
		double timesince = Math.min(dt,(L-dist)/(cv1-cv2));
		
		// move the boats backward in time
		b1.x -= b1.vx * timesince;
		b1.y -= b1.vy * timesince;
		b2.x -= b2.vx * timesince;
		b2.y -= b2.vy * timesince;
		*/
		
		// unit vector in the direction of the collision
		double ax = (b2.x-b1.x)/(L-overlap);
		double ay = (b2.y-b1.y)/(L-overlap);

		// velocity in the direction of the collision
		double cv1 = b1.vx*ax + b1.vy*ay;
		double cv2 = b2.vx*ax + b2.vy*ay;
		// velocity perpendicular to the collision
		double pv1 = -b1.vx*ay + b1.vy*ax;
		double pv2 = -b2.vx*ay + b2.vy*ax;

		// push the boats apart
		b1.x -= overlap*ax/2;
		b1.y -= overlap*ay/2;
		b2.x += overlap*ax/2;
		b2.y += overlap*ay/2;
		
		// keep the perpendicular components, swap the parallel
		if(cv1 - cv2 > 0){
			b1.vx = -pv1*ay + cv2*ax*restitution;
			b1.vy = pv1*ax + cv2*ay*restitution;
			b2.vx = -pv2*ay + cv1*ax*restitution;
			b2.vy = pv2*ax + cv1*ay*restitution;


			//b1.vx = (b2.vx*ax + b2.vy*ay)*ax - (-b1.vx*ay + b1.vy*ax)*ay;
			//b1.vy = (b2.vx*ax + b2.vy*ay)*ay + (-b1.vx*ay + b1.vy*ax)*ax;
			//b2.vx = (b1.vx*ax + b1.vy*ay)*ax - (-b2.vx*ay + b2.vy*ax)*ay;
			//b2.vy = (b1.vx*ax + b1.vy*ay)*ay + (-b2.vx*ay + b2.vy*ax)*ax;
		}
		
		/*
		// move time forward by timesince
		b1.x += b1.vx * timesince;
		b1.y += b1.vy * timesince;
		b2.x += b2.vx * timesince;
		b2.y += b2.vy * timesince;
		*/

		// recompute this distance
		distancesq[b1.index][b2.index] = (b1.x-b2.x)*(b1.x-b2.x) + (b1.y-b2.y)*(b1.y-b2.y);
		distancesq[b2.index][b1.index] = distancesq[b1.index][b2.index];
	}
	
	// sets up the given F matrix
	public double[][] setUpF(double[][] ff){
		// clear the vector
		for (int i=0; i<nx; i++){
			for (int j=0; j<ny; j++){
				ff[i][j] = 0.0;
			}
		}
		// for each boat, add a constant source at the relevant points
		int left, right, top, bottom;
		double distsq = 0.0;
		for(int i=0; i<N; i++){
			left = (int)((rk4Data[0][i][0] - L/2)/dx);
			right = (int)((rk4Data[0][i][0] + L/2)/dx)+1;
			bottom = (int)((rk4Data[0][i][1] - L/2)/dx);
			top = (int)((rk4Data[0][i][1] + L/2)/dx)+1;
			
			if(left < 0){ left += nx-1; }
			if(bottom < 0){ bottom = 0; }
			if(right > nx-1){ right -= nx-1; }
			if(top > ny-1){ top = ny-1; }
			
			// if the boat is not crossing the x boundary
			if(left < right){
				for(int xin=left; xin<=right; xin++){
					for(int yin=bottom; yin<=top; yin++){
						distsq = (rk4Data[0][i][0]-xin*dx)*(rk4Data[0][i][0]-xin*dx) + (rk4Data[0][i][1]-yin*dx)*(rk4Data[0][i][1]-yin*dx);
						if(distsq < rc*rc){
							ff[xin][yin] = a*(1.0-(c[xin][yin]/csat)*(c[xin][yin]/csat)*(c[xin][yin]/csat));
						}
						if(distsq < L*L/4){
							//ff[xin][yin] += c[xin][yin]*k;
						}
					}
				}
			}else{
				// if the boat position is on the left side of the screen
				if(rk4Data[0][i][0] < right*dx){
					for(int xin=0; xin<=right; xin++){
						for(int yin=bottom; yin<=top; yin++){
							distsq = (rk4Data[0][i][0]-xin*dx)*(rk4Data[0][i][0]-xin*dx) + (rk4Data[0][i][1]-yin*dx)*(rk4Data[0][i][1]-yin*dx);
							if(distsq < rc*rc){
								ff[xin][yin] = a*(1.0-(c[xin][yin]/csat)*(c[xin][yin]/csat)*(c[xin][yin]/csat));
							}
							if(distsq < L*L/4){
								//ff[xin][yin] += c[xin][yin]*k;
							}
						}
					}
					for(int xin=left; xin<=nx-1; xin++){
						for(int yin=bottom; yin<=top; yin++){
							distsq = (rk4Data[0][i][0]-(xin-(nx-1))*dx)*(rk4Data[0][i][0]-(xin-(nx-1))*dx) + (rk4Data[0][i][1]-yin*dx)*(rk4Data[0][i][1]-yin*dx);
							if(distsq < rc*rc){
								ff[xin][yin] = a*(1.0-(c[xin][yin]/csat)*(c[xin][yin]/csat)*(c[xin][yin]/csat));
							}
							if(distsq < L*L/4){
								//ff[xin][yin] += c[xin][yin]*k;
							}
						}
					}
				}else{
					for(int xin=0; xin<=right; xin++){
						for(int yin=bottom; yin<=top; yin++){
							distsq = (rk4Data[0][i][0]-(xin+(nx-1))*dx)*(rk4Data[0][i][0]-(xin+(nx-1))*dx) + (rk4Data[0][i][1]-yin*dx)*(rk4Data[0][i][1]-yin*dx);
							if(distsq < rc*rc){
								ff[xin][yin] = a*(1.0-(c[xin][yin]/csat)*(c[xin][yin]/csat)*(c[xin][yin]/csat));
							}
							if(distsq < L*L/4){
								//ff[xin][yin] += c[xin][yin]*k;
							}
						}
					}
					for(int xin=left; xin<=nx-1; xin++){
						for(int yin=bottom; yin<=top; yin++){
							distsq = (rk4Data[0][i][0]-xin*dx)*(rk4Data[0][i][0]-xin*dx) + (rk4Data[0][i][1]-yin*dx)*(rk4Data[0][i][1]-yin*dx);
							if(distsq < rc*rc){
								ff[xin][yin] = a*(1.0-(c[xin][yin]/csat)*(c[xin][yin]/csat)*(c[xin][yin]/csat));
							}
							if(distsq < L*L/4){
								//ff[xin][yin] += c[xin][yin]*k;
							}
						}
					}
				}
				for(int xin=0; xin<=right; xin++){
					for(int yin=bottom; yin<=top; yin++){
						distsq = (rk4Data[0][i][0]-xin*dx)*(rk4Data[0][i][0]-xin*dx) + (rk4Data[0][i][1]-yin*dx)*(rk4Data[0][i][1]-yin*dx);
						if(distsq < rc*rc){
							ff[xin][yin] = a*(1.0-(c[xin][yin]/csat)*(c[xin][yin]/csat)*(c[xin][yin]/csat));
						}
						if(distsq < L*L/4){
							//ff[xin][yin] += c[xin][yin]*k;
						}
					}
				}
				for(int xin=left; xin<=nx-1; xin++){
					for(int yin=bottom; yin<=top; yin++){
						distsq = (rk4Data[0][i][0]-xin*dx)*(rk4Data[0][i][0]-xin*dx) + (rk4Data[0][i][1]-yin*dx)*(rk4Data[0][i][1]-yin*dx);
						if(distsq < rc*rc){
							ff[xin][yin] = a*(1.0-(c[xin][yin]/csat)*(c[xin][yin]/csat)*(c[xin][yin]/csat));
						}
						if(distsq < L*L/4){
							//ff[xin][yin] += c[xin][yin]*k;
						}
					}
				}
			}
		}
		
		return ff;
	}

	public void solveSystem(int direction, int index){
		// solves the tridiagonal part
		// for the boundary: do the first and last points separately
		// diag -> kdt/4 , upper=lower -> 0
		// direction: 1=x, -1=y
		if(direction == 1){
			// periodic x direction requires (sigh)
			solveXDirection(index);
			/*
			// the x direction, first index
			for (int i=1; i<nx-1; i++){
				//tmpm = lower/diag;
				//tmpdiag = diag - tmpm*upper;
				rhs[i] = rhs[i] - tmpm*rhs[i-1];
			}
			rhs[nx-1] = rhs[nx-1] - 2*tmpm*rhs[nx-2]; // this is for the boundary
			
			tmpc[nx-1][index] = rhs[nx-1]/(diag - 2*tmpm*upper); // this is for the boundary
			for (int i=nx-2; i>1; i--){
				tmpc[i][index] = (rhs[i]-upper*tmpc[i+1][index])/tmpdiag;
			}
			tmpc[1][index] = (rhs[1]-upper*tmpc[2][index])/(diag - tmpm*2*upper); // this is for the boundary
			tmpc[0][index] = (rhs[0]-2*upper*tmpc[1][index])/diag; // this is for the boundary
			*/
		}else{
			// the y direction, second index
			for (int i=1; i<ny-1; i++){
				//tmpm = lower/diag;
				//tmpdiag = diag - tmpm*upper;
				rhs[i] = rhs[i] - tmpm*rhs[i-1];
			}
			rhs[ny-1] = rhs[ny-1] - 2*tmpm*rhs[ny-2]; // this is for the boundary

			tmpc[index][ny-1] = rhs[ny-1]/(diag - 2*tmpm*upper); // this is for the boundary
			for (int i=ny-2; i>1; i--){
				tmpc[index][i] = (rhs[i]-upper*tmpc[index][i+1])/tmpdiag;
			}
			tmpc[index][1] = (rhs[1]-upper*tmpc[index][2])/(diag - tmpm*2*upper); // this is for the boundary
			tmpc[index][0] = (rhs[0]-2*upper*tmpc[index][1])/diag; // this is for the boundary
		}
	}
	
	public void solveXDirection(int index){
		// uses sherman morrison method
		// solve Bc = d by By=d, Bq=u, c=y-(vy)/(1+vq)q
		// u = [-1..lower]
		// v = [1..-upper]    (maybe change some signs if it doesn't work)
		double id = 0.0;
		// solve for y
		mc[0] = upper/(diag+1);
		md[0] = rhs[0]/(diag+1);
		for (int i = 1; i < nx-1; i++){
			id = 1/(diag - mc[i-1]*lower);
			mc[i] = upper*id;
			md[i] = (rhs[i] - md[i-1]*lower)*id;
        }
		id = 1/(diag + upper*lower - mc[nx-2]*lower);
		mc[nx-1] = upper*id;
		md[nx-1] = (rhs[nx-1] - md[nx-2]*lower)*id;
		
        y[nx-1] = md[nx-1];
        for (int i=nx-2; i>=0; i--){
			y[i]=md[i]-mc[i]*y[i+1];
		}
		
		//solve for q
		mc[0] = upper/(diag+1);
		md[0] = -1.0/(diag+1);
		for (int i = 1; i < nx-1; i++){
			id = 1/(diag - mc[i-1]*lower);
			mc[i] = upper*id;
			md[i] = -md[i-1]*lower*id;
        }
		id = 1/(diag + upper*lower - mc[nx-2]*lower);
		mc[nx-1] = upper*id;
		md[nx-1] = (lower - md[nx-2]*lower)*id;
		
        q[nx-1] = md[nx-1];
        for (int i=nx-2; i>=0; i--){
			q[i]=md[i]-mc[i]*q[i+1];
		}
		
		//find the products
		double vy = y[0]-y[nx-1]*upper;
		double vq = q[0]-q[nx-1]*upper;
		
		//finally solve for c
		for (int i = 0; i < nx; i++){
			tmpc[i][index] = y[i] - (vy/(1+vq))*q[i];
        }
	}
	
	///////////////////////////////////////////////////////////////////////////
	// this is the HYBRID time integrator 
	// It does: 1. RK4 on motion eq. for a half time step (for fnh)
	//			2. RK4 again for the second half step (for fnp)
	//			3. ADI full step
	// 
	// it does 1 time step per call
	///////////////////////////////////////////////////////////////////////////
	public void timeStep(){
		
		
		// perform a half time step
		stepRK4(dt/2, false);
		
		// set up fnh
		fnh = setUpF(fnh);
		
		// perform the full step
		stepRK4(dt, false);
		
		// set up fnp
		fnp = setUpF(fnp);
		
		// do the ADI step
		fullStepADI();
		
		// do RK4 again with the updated data
		stepRK4(dt, true);
		
		// set up fn
		fn = setUpF(fn);
		
		// compute stats
		if(computeStats && currentT%4 == 0){
			computeStatistics((int)(currentT/4));
		}

		//advance time
		currentT++;
	}
	
	// performs one RK4 half time step on the motion equations
	public void stepRK4(double step, boolean complete){
		Buoy2D bn;
		// perform a time step
		//prepare temp0 and step1
		for(int i=0; i<N; i++){
			bn = b.get(i);
			rk4Data[0][i][0] = bn.x;
			rk4Data[0][i][1] = bn.y;
			rk4Data[0][i][2] = bn.vx;
			rk4Data[0][i][3] = bn.vy;
			
			rk4Data[1][i][0] = bn.x;
			rk4Data[1][i][1] = bn.y;
			rk4Data[1][i][2] = bn.vx;
			rk4Data[1][i][3] = bn.vy;
		}
		rk4Data[1] = motionEquation(rk4Data[1]);	//first
		for(int i=0; i<N; i++){
			rk4Data[1][i][0] = rk4Data[1][i][0]*step;
			rk4Data[1][i][1] = rk4Data[1][i][1]*step;
			rk4Data[1][i][2] = rk4Data[1][i][2]*step;
			rk4Data[1][i][3] = rk4Data[1][i][3]*step;
		}
		
		for(int i=0; i<N; i++){
			rk4Data[2][i][0] = rk4Data[0][i][0] + rk4Data[1][i][0]/2;
			rk4Data[2][i][1] = rk4Data[0][i][1] + rk4Data[1][i][1]/2;
			rk4Data[2][i][2] = rk4Data[0][i][2] + rk4Data[1][i][2]/2;
			rk4Data[2][i][3] = rk4Data[0][i][3] + rk4Data[1][i][3]/2;
		}
		rk4Data[2] = motionEquation(rk4Data[2]);	//second
		for(int i=0; i<N; i++){
			rk4Data[2][i][0] = rk4Data[2][i][0]*step;
			rk4Data[2][i][1] = rk4Data[2][i][1]*step;
			rk4Data[2][i][2] = rk4Data[2][i][2]*step;
			rk4Data[2][i][3] = rk4Data[2][i][3]*step;
		}
		
		for(int i=0; i<N; i++){
			rk4Data[3][i][0] = rk4Data[0][i][0] + rk4Data[2][i][0]/2;
			rk4Data[3][i][1] = rk4Data[0][i][1] + rk4Data[2][i][1]/2;
			rk4Data[3][i][2] = rk4Data[0][i][2] + rk4Data[2][i][2]/2;
			rk4Data[3][i][3] = rk4Data[0][i][3] + rk4Data[2][i][3]/2;
		}
		rk4Data[3] = motionEquation(rk4Data[3]);	//third
		for(int i=0; i<N; i++){
			rk4Data[3][i][0] = rk4Data[3][i][0]*step;
			rk4Data[3][i][1] = rk4Data[3][i][1]*step;
			rk4Data[3][i][2] = rk4Data[3][i][2]*step;
			rk4Data[3][i][3] = rk4Data[3][i][3]*step;
		}
		
		for(int i=0; i<N; i++){
			rk4Data[4][i][0] = rk4Data[0][i][0] + rk4Data[3][i][0];
			rk4Data[4][i][1] = rk4Data[0][i][1] + rk4Data[3][i][1];
			rk4Data[4][i][2] = rk4Data[0][i][2] + rk4Data[3][i][2];
			rk4Data[4][i][3] = rk4Data[0][i][3] + rk4Data[3][i][3];
		}
		rk4Data[4] = motionEquation(rk4Data[4]);	//fourth
		for(int i=0; i<N; i++){
			rk4Data[4][i][0] = rk4Data[4][i][0]*step;
			rk4Data[4][i][1] = rk4Data[4][i][1]*step;
			rk4Data[4][i][2] = rk4Data[4][i][2]*step;
			rk4Data[4][i][3] = rk4Data[4][i][3]*step;
		}
		
		// store the step in temp0
		for(int i=0; i<N; i++){
			rk4Data[0][i][0] += (rk4Data[1][i][0] + 2*rk4Data[2][i][0] + 2*rk4Data[3][i][0] + rk4Data[4][i][0])/6;
			rk4Data[0][i][1] += (rk4Data[1][i][1] + 2*rk4Data[2][i][1] + 2*rk4Data[3][i][1] + rk4Data[4][i][1])/6;
			rk4Data[0][i][2] += (rk4Data[1][i][2] + 2*rk4Data[2][i][2] + 2*rk4Data[3][i][2] + rk4Data[4][i][2])/6;
			rk4Data[0][i][3] += (rk4Data[1][i][3] + 2*rk4Data[2][i][3] + 2*rk4Data[3][i][3] + rk4Data[4][i][3])/6;
		}
		
		if(complete){
			// if this is the final step to be stored
			for(int i=0; i<N; i++){
				bn = b.get(i);
				bn.x = rk4Data[0][i][0];
				bn.y = rk4Data[0][i][1];
				bn.vx = rk4Data[0][i][2];
				bn.vy = rk4Data[0][i][3];
				
				if(withNoise){
					bn.vx += noiseStrength*dt*noiseMaker.nextGaussian();
					bn.vy += noiseStrength*dt*noiseMaker.nextGaussian();
				}
			}
			
			computeDistances();
			doCollisions();
			
		}
		
	}
	
	// performs one full ADI step
	public void fullStepADI(){
		////////////////////////////////////////////////////////
		// do one step of the ADI method (x direction then y)
		// do it one column(y value) at a time, then one row(x value) at a time
		// the boundary condition is supposed to be grad(c)*normal = 0
		// who knows if it will really come out that way... (you aren't supposed to say that! What are you doing!)
		////////////////////////////////////////////////////////
		//////////     Ew=d     ///////////
		// lam = D*dt/2/dx^2
		// wy[k] = c[k][y]
		//
		// d = lam*(wy(k-1) + wy(k+1)) + (1-2lam-k*dt/4)*wy(k) + by(k) + a*dt/4(Fn+Fnh)
		// B = [Ex+(k*dt/4)*I] (already made)
		//
		// d = lam*(wx(k-1) + wx(k+1)) + (1-2lam-k*dt/4)*wx(k) + bx(k) + a*dt/4(Fnh+Fnp)
		// B = [Ey+(k*dt/4)*I] (already made)
		//
		////////////////////////////////////////////////////////
		
		// do it one yind at a time ///////////////////////////////////////////
		for(int yind=1; yind<ny-1; yind++){
			// set up d 
			for (int i=0; i<nx; i++){
				rhs[i] = lambda*(c[i][yind-1]+c[i][yind+1]) + usefulnumber*c[i][yind] + adt*(fn[i][yind]+fnh[i][yind]);
			}
			//B is already set, so solve the system for c
			solveSystem(1, yind);
		}
		// and for the boundaries
		// no flux boundary
		for (int i=0; i<nx; i++){
			rhs[i] = lambda*(2*c[i][1]) + usefulnumber*c[i][0] + adt*(fn[i][0]+fnh[i][0]);
		}
		solveSystem(1, 0);
		for (int i=0; i<nx; i++){
			rhs[i] = lambda*(2*c[i][ny-2]) + usefulnumber*c[i][ny-1] + adt*(fn[i][ny-1]+fnh[i][ny-1]);
		}
		solveSystem(1, ny-1);
		
		// transfer result from tmpc to c
		for (int i=0; i<nx; i++){
			for (int j=0; j<ny; j++){
				c[i][j] = tmpc[i][j];
			}
		}
		
		
		// do it one xind at a time ///////////////////////////////////////////
		for(int xind=1; xind<nx-1; xind++){
			// set up d 
			for (int i=0; i<ny; i++){
				rhs[i] = lambda*(c[xind-1][i]+c[xind+1][i]) + usefulnumber*c[xind][i] + adt*(fnh[xind][i]+fnp[xind][i]);
			}
			//B is already set, so solve the system for c
			solveSystem(-1, xind);
		}
		// and for the boundaries
		// Periodic boundary
		for (int i=0; i<ny; i++){
			rhs[i] = lambda*(c[1][i] + c[nx-1][i]) + usefulnumber*c[0][i] + adt*(fnh[0][i]+fnp[0][i]);
		}
		solveSystem(-1, 0);
		for (int i=0; i<ny; i++){
			rhs[i] = lambda*(c[nx-2][i] + c[0][i]) + usefulnumber*c[nx-1][i] + adt*(fnh[nx-1][i]+fnp[nx-1][i]);
		}
		solveSystem(-1, nx-1);
		
		// transfer result from tmpc to c
		for (int i=0; i<nx; i++){
			for (int j=0; j<ny; j++){
				c[i][j] = tmpc[i][j];
				// avoid negative c if needed
				if(c[i][j] < 0){
					//c[i][j] = 0.0;
				}
			}
		}
		
	}
	
	///////////////////////////////////////////////////////////////////////////
	// Thus ends the integrator
	///////////////////////////////////////////////////////////////////////////
	
	// computes various statistics
	
	public void computeStatistics(int index){
		/*
		if(index >= meanDisp.length){
			index = index%meanDisp.length;
		}
		// mean displacement
		double dx, dy;
		double tmpdispx = 0;
		double tmpdispy = 0;
		double tmpr2 = 0;
		for(int i=0; i<N; i++){
			dx = b.get(i).x - initialx[i];
			dy = b.get(i).y - initialy[i];
			tmpdispx += dx;
			tmpdispy += dy;
			tmpr2 += dx*dx+dy*dy;
		}
		meanDisp[index] = Math.sqrt(tmpdispx*tmpdispx/N/N + tmpdispy*tmpdispy/N/N);
		meansqDisp[index] = tmpr2/N;
		
		// mean square displacement
		
		// average minimum distance
		double tmpmindist = 0;
		currentMinDist = 0;
		for(int i=0; i<N; i++){
			tmpmindist = R;
			for(int j=0; j<N; j++){
				if(distancesq[i][j] < tmpmindist*tmpmindist && i!=j){
					tmpmindist = Math.sqrt(distancesq[i][j]);
				}
			}
			tmpmindist -= L;
			currentMinDist += tmpmindist;
		}
		currentMinDist = currentMinDist/N;
		aveminDist[index] = currentMinDist;
		
		// average speed
		avev[index] = 0.0;
		for(int i=0; i<N; i++){
			avev[index] += Math.sqrt(b.get(i).vx*b.get(i).vx + b.get(i).vy*b.get(i).vy);
		}
		avev[index] = avev[index]/N;
		 */
	}
	 

	/**
	 * adds one boat to the system in some available space
	 */
	public int addBoat(){
		int placed = -1;
		// try to stick it in the system
		if(N == 0){
			// first boat goes in the center;
			Buoy2D newb = new Buoy2D(0);
			newb.x = R/2;
			newb.y = W/2;
			newb.vx = 0.0;
			newb.vy = 0.0;
			
			b.add(0, newb);
			N++;
			placed = 0;
		}else {
			// put it in a random spot
			Buoy2D newb = new Buoy2D(N);
			Random rand = new Random();
			newb.x = rand.nextDouble()*(R-L-L) + L;
			newb.y = rand.nextDouble()*(W-L-L) + L;
			newb.vx = 0.0;
			newb.vy = 0.0;
			b.add(N, newb);
			N++;
			placed = N;
		}
		if(placed >= 0){
			// take care of some allocation details
			// keeping the old data.
			double[][][] temprk4 = new double[5][N][4];
			for(int p=0; p<5; p++){
				for(int i=0; i<N-1; i++){
					for(int j=0; j<4; j++){
						temprk4[p][i][j] = rk4Data[p][i][j];
					}
				}
			}
			rk4Data = temprk4;

			netx = new double[N];
			nety = new double[N];
			
			distancesq = new double[N][N];

			initialx = new double[N];
			initialy = new double[N];
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
		
		// notice I didn't deallocate...
		
		return N;
	}

	public void homogenize(){
		if(N>0){
			// sorry, this doesn't do anything yet
		}
	}

	public void randomize(){
		if(N>0){
			Random rand = new Random();
			for(int i=0; i<N; i++){
				b.get(i).x = rand.nextDouble()*(R-L-L) + L;
				b.get(i).y = rand.nextDouble()*(W-L-L) + L;
			}
			computeDistances();
			doCollisions();
		}
	}
	
	public int getN(){
		return N;
	}

	public ArrayList<Buoy2D> getBoats(){
		return b;
	}

	public double[][] getc(){
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
		nt = (long)(T*stepspersec);
	}

	public void setN(int newN){
		reset();
		for(int i=0; i<newN; i++){
			addBoat();
		}
	}
	
	public void setRandW(double newR, double newW){
		nx = (int)(stepsperrc*newR/rc);
		ny = (int)(stepsperrc*newW/rc);
		R = newR;
		W = newW;
		dx = R/(nx-1);
		xg = new double[nx];
		yg = new double[ny];
		for(int i=0; i<nx; i++){
			xg[i] = dx*i;
		}
		for(int i=0; i<ny; i++){
			yg[i] = dx*i;
		}
		
		c = new double[nx][ny];

		netx = new double[N];
		nety = new double[N];
		
		rk4Data = new double[5][N][4];
		
		tmpc = new double[nx][ny];
		rhs = new double[nx];
		fn = new double[nx][ny];
		fnh = new double[nx][ny];
		fnp = new double[nx][ny];
		
		mc = new double[nx];
		md = new double[nx];
		y = new double[nx];
		q = new double[nx];
		
		lambda = D*dt/2/dx/dx;
		usefulnumber = 1-2*lambda - dt*k/4;
		adt = dt/4;
		diag = 1.0+2*lambda + dt*k/4;
		upper = -lambda;
		lower = -lambda;
		tmpm = lower/diag;
		tmpdiag = diag - tmpm*upper;
		
		//initialize c
		for(int i=0; i<nx; i++){
			for(int j=0; j<ny; j++){
				c[i][j] = 0;
			}
		}
	}

	public double[] getParameters(){
		double[] params = new double[13];
		params[0] = R;
		params[1] = W;
		params[2] = m;
		params[3] = vis;
		params[4] = L;
		params[5] = rc;
		params[6] = D;
		params[7] = k;
		params[8] = a;
		params[9] =  beta;
		params[10] = csat;
		params[11] = maxdgamma;
		params[12] = restitution;
		return params;
	}

	public void setParameters(double[] params){
		if((R != params[0]) || (W != params[1])){
			setRandW(params[0], params[1]);
		}
		m = params[2];
		vis = params[3];
		L = params[4];
		rc = params[5];
		D = params[6];
		k = params[7];
		a = params[8];
		beta = params[9];
		csat = params[10];
		maxdgamma = params[11];
		restitution = params[12];
		
		// set some important stuff
		lambda = D*dt/2/dx/dx;
		usefulnumber = 1-2*lambda - dt*k/4;
		adt = dt/4;
		diag = 1.0+2*lambda + dt*k/4;
		upper = -lambda;
		lower = -lambda;
		tmpm = lower/diag;
		tmpdiag = diag - tmpm*upper;
	}

	public double[] getDisc(){
		double[] disc = new double[8];
		disc[0] = T;
		disc[1] = dt;
		disc[2] = (double)nt;
		disc[3] = R;
		disc[4] = W;
		disc[5] = dx;
		disc[6] = (double)nx;
		disc[7] = (double)ny;
		return disc;
	}
	
	public void setDisc(int nxperrc, int ntpersec){
		stepsperrc = nxperrc;
		nx = (int)(stepsperrc*R/rc);
		ny = (int)(stepsperrc*W/rc);
		dx = R/(nx-1);
		xg = new double[nx];
		yg = new double[ny];
		for(int i=0; i<nx; i++){
			xg[i] = dx*i;
		}
		for(int i=0; i<ny; i++){
			yg[i] = dx*i;
		}
		
		stepspersec = ntpersec;
		currentT = 0;
		nt = (long)(T*stepspersec);
		dt = T/(nt-1);
		
		c = new double[nx][ny];
		
		tmpc = new double[nx][ny];
		rhs = new double[nx];
		fn = new double[nx][ny];
		fnh = new double[nx][ny];
		fnp = new double[nx][ny];
		
		mc = new double[nx];
		md = new double[nx];
		y = new double[nx];
		q = new double[nx];
		
		lambda = D*dt/2/dx/dx;
		usefulnumber = 1-2*lambda - dt*k/4;
		adt = dt/4;
		diag = 1.0+2*lambda + dt*k/4;
		upper = -lambda;
		lower = -lambda;
		tmpm = lower/diag;
		tmpdiag = diag - tmpm*upper;
		
		//initialize c
		for(int i=0; i<nx; i++){
			for(int j=0; j<ny; j++){
				c[i][j] = 0;
			}
		}
	}

	public void reset(){
		currentT = 0;
		for(int i=N; i>0; i--){
			removeBoat();
		}
		for(int i=0; i<nx; i++){
			for(int j=0; j<ny; j++){
				c[i][j] = 0;
			}
		}
	}

	public void startRecording(){
		computeStats = true;
		currentT = 0;
		for(int i=0; i<N; i++){
			initialx[i] = b.get(i).x;
			initialy[i] = b.get(i).y;
		}
		meanDisp = new double[(int)(nt/4)];
		meansqDisp = new double[(int)(nt/4)];
		aveminDist = new double[(int)(nt/4)];
		currentMinDist = 0;
		avev = new double[(int)(nt/4)];
	}

	public void stopRecording(){
		computeStats = false;
	}
	
	public void changeBackground(double backc){
		background = backc;
	}
	
	public void kickABoat(int ind, double kick){
		b.get(ind).vx += kick;
	}
	
	public void addNoise(boolean addit, double amount){
		withNoise = addit;
		noiseStrength = amount;
		noiseMaker = new Random();
	}
}
