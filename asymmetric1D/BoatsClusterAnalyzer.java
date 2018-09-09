/*
 * stores some data and does some math with it
 */
public class BoatsClusterAnalyzer {
	double[][] xdata, frontdata;
	double[][] beta;
	int N, points, position;
	double dt, D;

	public BoatsClusterAnalyzer(int nN, int npoints, double ndt, double nD){
		N = nN;
		points = npoints;
		dt = ndt;
		D = nD;
		
		position = 0;
		xdata = new double[points][N];
		frontdata = new double[points][N];
	}
	
	public void addData(double[] x, double[] fr){
		if(position >= points){
			return;
		}
		for(int i=0; i<N; i++){
			xdata[position][i] = x[i];
			frontdata[position][i] = fr[i];
		}
		
		position++;
	}
	
	public void check(){
		for(int i=0; i< position; i++){
			System.out.println("x= "+String.valueOf(xdata[i][0]));
		}
	}
	
	// beta is the normalized instantaneous difference between front space and back space
	// computed over (timesAround) cycles around the circle
	// D is the nomalization constant
	public void computeBeta(){
		
		int cycle = points-2;
		int timesAround = 0;
		// count the cycles and find the end of an integer number of them
		for(int i=1; i<points-1; i++){
			if(xdata[i][0] < xdata[0][0] && xdata[i+1][0] >= xdata[0][0]){
				timesAround++;
				cycle = i;
			}
		}
		
		// compute beta for timesAround cycles
		if(timesAround < 1){
			System.out.println("didn't record long enough");
		}
		
		beta = new double[cycle][N];
		for(int i=0; i<cycle; i++){
			for(int j=1; j<N; j++){
				beta[i][j] = (frontdata[i][j] - frontdata[i][j-1])*D;
			}
			beta[i][0] = (frontdata[i][0] - frontdata[i][N-1])*D;
		}
	}
	
	// computes the time average of beta for each boat
	public double[] computeAveBeta(){
		//compute the average beta for each boat
		computeBeta();
		double[] avebeta = new double[N];
		
		for(int i=0; i<N; i++){
			avebeta[i] = 0.0;
			for(int j=0; j<beta.length; j++){
				avebeta[i] += beta[j][i];
			}
			avebeta[i] = avebeta[i]/beta.length;
		}
		
		return avebeta;
	}
	
	// returns the max average value of beta
	public double computeMaxBeta(){
		
		double[] avebeta = computeAveBeta();
		
		double maxbeta = -1;
		for(int i=0; i<N; i++){
			maxbeta = Math.max(maxbeta,avebeta[i]);
		}
		
		return maxbeta;
	}
	
	// plots the max average beta vs time
	public void plotBetaOfT(){
		if(beta == null){
			return;
		}
		int timeskip = 100;
		int tpts = beta.length/timeskip;
		// the time data
		double[] tdata = new double[tpts];
		int k = 0;
		for(int i=0; i<beta.length-timeskip; i+=timeskip){
			tdata[k] = i*dt;
			k++;
		}
		// the max beta data
		double[] bdata = new double[tpts];
		int maxboat;
		double[] avebeta;
		for(int i=0; i<tpts; i++){
			avebeta = computeAveBeta(tdata[i]);
			maxboat = 0;
			for(int j=1; j<N; j++){
				if(avebeta[j] > avebeta[j-1]){
					maxboat = j;
				}
			}
			bdata[i] = avebeta[maxboat];
		}
		
		Plot2D betaplot = new Plot2D(700, 500);
		betaplot.showAxes(true);
		betaplot.setLabels("time", "frontspace-backspace");
		betaplot.boxed(true);
		betaplot.lines(true);
		betaplot.addX(tdata);
		betaplot.addY(bdata);
		betaplot.setXLimits(0.0, tdata[tdata.length-1]);
		
		betaplot.paint();
		
		PlotWindow pwind = new PlotWindow(betaplot);
	}
	
	// computes the PARTIAL time average of beta for each boat
	public double[] computeAveBeta(double maxt){
		//compute the average beta for each boat
		computeBeta();
		double[] avebeta = new double[N];
		
		for(int i=0; i<N; i++){
			avebeta[i] = 0.0;
			for(int j=0; j<beta.length; j++){
				if(j*dt > maxt){ break; }
				avebeta[i] += beta[j][i];
			}
			avebeta[i] = avebeta[i]/beta.length;
		}
		
		return avebeta;
	}
}
