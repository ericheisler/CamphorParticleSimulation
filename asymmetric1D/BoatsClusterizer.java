import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.util.*;
import java.io.*;
import javax.imageio.ImageIO;

/**
 *
 * @author Bilbo
 */
public class BoatsClusterizer extends JPanel implements ActionListener{

	Vector<Boat> boats;
	double[] x, c, y, q, aF, mc, md, force;
	double netForce;
	double diag, upper, lower;
	double R, m, vis, l, L, rc, D, k, a, b, v;
	int nx;
	double dx;
	boolean vset;

	int width, height;
	JFrame frame;
	Plot2D plot;
	PPanel parPan;
	JPanel controlPanel;
	JButton calcButton, adButton, rmButton, saveButton, setButton;
	JButton setvButton, shakeButton, testButton, seriesButton, nseriesButton, infButton;


	public BoatsClusterizer(double[] params){
		super();
		setR(params[0]);
		m = params[1];
		vis = params[2];
		l = params[3];
		L = params[4];
		rc = params[5];
		D = params[6];
		k = params[7];
		a = params[8];
		b = params[9];

		v = 30;
		netForce = 0.0;
		vset = false;
		force = new double[0];

		boats = new Vector<Boat>();

		diag = 2*(D/dx/dx)+k;
		upper = -D/dx/dx - v/dx;
		lower = -D/dx/dx + v/dx;

		width = 800;
		height = 450;
		setPreferredSize(new Dimension(width, height));
		plot = new Plot2D(width, height-50);
		plot.showAxes(true);
		plot.setLabels("x", "c");
		plot.boxed(true);
		plot.lines(true);
		plot.paint();

		parPan = new PPanel(params);

		//create some buttons
		calcButton = new JButton("calculate");
		calcButton.addActionListener(this);
		adButton = new JButton("add boat");
		adButton.addActionListener(this);
		rmButton = new JButton("remove boat");
		rmButton.addActionListener(this);
		saveButton = new JButton("save plot");
		saveButton.addActionListener(this);
		setButton = new JButton("set parameters");
		setButton.addActionListener(this);
		setvButton = new JButton("set velocity");
		setvButton.addActionListener(this);
		shakeButton = new JButton("shake it up");
		shakeButton.addActionListener(this);
		testButton = new JButton("run test boat");
		testButton.addActionListener(this);
		seriesButton = new JButton("vary pars");
		seriesButton.addActionListener(this);
		nseriesButton = new JButton("vary N");
		nseriesButton.addActionListener(this);
		infButton = new JButton("to infinity!!");
		infButton.addActionListener(this);

		controlPanel = new JPanel();
		controlPanel.setLayout(new BoxLayout(controlPanel, BoxLayout.LINE_AXIS));
		controlPanel.add(calcButton);
		controlPanel.add(adButton);
		controlPanel.add(rmButton);
		controlPanel.add(saveButton);
		controlPanel.add(setButton);
		controlPanel.add(setvButton);
		controlPanel.add(shakeButton);
		controlPanel.add(testButton);
		controlPanel.add(seriesButton);
		controlPanel.add(nseriesButton);
		controlPanel.add(infButton);

		frame = new JFrame("boats simulation");
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

		Container pane = frame.getContentPane();
		pane.setLayout(new BorderLayout());
		pane.add(this, BorderLayout.CENTER);
		pane.add(controlPanel, BorderLayout.SOUTH);
		pane.add(parPan, BorderLayout.EAST);

		pane.validate();
		//Display the window.
        frame.pack();
        frame.setVisible(true);

	}

	// for the actionlistener
	public void actionPerformed(ActionEvent e){
		if (e.getSource() == calcButton){
			calc();
		} else if (e.getSource() == adButton){
			addBoat();
		} else if (e.getSource() == rmButton){
			removeBoat();
		} else if (e.getSource() == saveButton){
			writeImage();
		} else if (e.getSource() == setButton){
			parPan.newParameters();
			diag = 2*(D/dx/dx)+k;
			upper = -D/dx/dx + v/dx;
			lower = -D/dx/dx - v/dx;
		} else if (e.getSource() == setvButton){
			if(!vset){
				vset = true;
				setvButton.setText("unset velocity");
				v = Double.valueOf((String)JOptionPane.showInputDialog(frame, "Input constant velocity"));
				diag = 2*(D/dx/dx)+k;
				upper = -D/dx/dx - v/dx;
				lower = -D/dx/dx + v/dx;
			}else{
				vset = false;
				setvButton.setText("set velocity");
			}
		} else if (e.getSource() == shakeButton){
			Random rand = new Random();
			for(Boat bn: boats){
				double dis = rand.nextDouble()*0.1;
				for(int i=bn.index; i<boats.size(); i++){
					boats.get(i).x -= dis;
				}
			}
		} else if (e.getSource() == testButton){
			testBoat();
		} else if (e.getSource() == seriesButton){
			testSeries();
		} else if (e.getSource() == nseriesButton){
			testNSeries();
		} else if (e.getSource() == infButton){
			toInfinity();
		}
	}

	// ye olde painting method
	public void paintComponent(Graphics g) {
		super.paintComponent(g);

		plot.clearData();
		plot.addX(x);
		plot.addY(c, plot.DOT, Color.RED);
		for(int j=0; j<force.length; j++){
			plot.addPoint(boats.get(j).x, force[j], plot.RING, Color.BLUE); 
		}
		plot.paint();
		g.drawImage(plot, 0, 0, this);

		g.setColor(Color.BLACK);
		g.drawString("v = "+String.format("%.6f", v), 25, 20);
		g.drawString("f = "+String.format("%.6f", netForce), 25, 35);
		
		//draw the boats
		g.setColor(Color.BLACK);
		int offset = 20;
		int csize = (int)(2*rc/R*(width-offset));
		int bsize = (int)(L/R*(width-offset));
		for(Boat bn: boats){
			int bpos = (int) (bn.x/R * (width-offset));
			int cpos = (int) ((bn.x-rc)/R * (width-offset));
			g.fillOval(bpos+offset, height-50+(int)(50/2 - bsize/2), bsize, bsize);
			g.drawOval(cpos+offset, height-50+(int)(50/2 - csize/2), csize, csize);
			if(bn.index == 0){
				g.setColor(Color.RED);
				g.fillOval(bpos+offset, height-50+(int)(50/2 - bsize/2), bsize, bsize);
				g.setColor(Color.BLACK);
			}
		}
	}
	
	public void calc(){
		if(boats.isEmpty()){return;}
		setaF();
		solvec();
		solveForce();
		repaint();
		
		//now iteratively relax the system
		boolean equilibrium = false;
		int iter = 0;
		int maxiters = 100;
		Random rand = new Random();
		
		while(!equilibrium && iter < maxiters){
			equilibrium = false;
			iter++;
			
			//if the boat in front feels larger force, increase separation
			//if it feels smaller force, decrease (until touching)
			double groupForce = force[0];
			int groupsize = 1;
			for(Boat bn: boats){
				if(bn.index == 0){
					//do nothing 
				}else{
					if(bn.x >= boats.get(bn.index-1).x-L-rc && groupForce/groupsize < force[bn.index]){
						groupForce += force[bn.index];
						groupsize++;
					}
					else{
						break;
					}
				}
			}
			double fmin = 1000;
			double fmax = -1000;
			for(int j=0; j<force.length; j++){
				fmin = Math.min(fmin, force[j]);
				fmax = Math.max(fmax, force[j]);
			}
			force[0] = groupForce/groupsize;
			for(Boat bn: boats){
				if(bn.index == 0){
					//do nothing because I don't want to move it
				}else{
					if(bn.x >= boats.get(bn.index-1).x-L-rc && groupForce/groupsize < force[bn.index]){
						force[bn.index] = groupForce/groupsize;
					}
					bn.x += (force[bn.index]-force[bn.index-1])*0.2/(fmax-fmin);
					bn.x = Math.min(bn.x, boats.get(bn.index-1).x-L-rc);
				}
			}
			
			// find the average force
			netForce = 0.0;
			for(int j=0; j<force.length; j++){
				netForce += force[j]; 
			}
			netForce = netForce/force.length;
			fmin = 1000;
			fmax = -1000;
			for(int j=0; j<force.length; j++){
				fmin = Math.min(fmin, force[j]);
				fmax = Math.max(fmax, force[j]);
			}
			//adjust the system velocity
			if(!vset){
				v += netForce*10;
				upper = -D/dx/dx - v/dx;
				lower = -D/dx/dx + v/dx;
			}
			
			int rem = 0;
			for(Boat bn: boats){
				if(bn.x <= 0){
					rem++;
				}
			}
			for(int j=0; j<rem; j++){
				removeBoat();
			}
			setaF();
			solvec();
			solveForce();
			
			paintImmediately(0, 0, width, height);
		}
	}
	
	////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////
	
	
	public void testBoat(){
		if(boats.size() < 1){
			return;
		}
		addBoat();
		int testind = boats.size()-1;
		double testdx = 0.1;
		int testpoints = 100;
		double[] testforce = new double[testpoints];
		double[] spacing = new double[testpoints];
		
		//move the test boat and check force
		boats.get(testind).x = boats.get(testind-1).x - L - rc;
		for(int i=0; i< testpoints; i++){
			boats.get(testind).x -= testdx;
			
			setaF();
			solvec();
			solveForce();
			
			testforce[i] = force[testind] - force[testind-1];
			spacing[i] = testdx*(i+1);
		}
		
		//do something with the data
		Plot2D testplot = new Plot2D(width, height-50);
		testplot.showAxes(true);
		testplot.setLabels("dist. from next boat", "force on test boat");
		testplot.setTitle("D="+String.format("%.2f", D)+" k="+String.format("%.4f", k));
		testplot.boxed(true);
		testplot.lines(true);
		testplot.addX(spacing);
		testplot.addY(testforce);
		testplot.setXLimits(0.0, spacing[testpoints-1]+testdx);
		testplot.setYLimits(-0.25, 0.25);
		
		testplot.paint();
		
		PlotWindow pwind = new PlotWindow(testplot);
		
		removeBoat();
	}
	
	public void testSeries(){
		if(boats.size() < 1){
			return;
		}
		Plot2D testplot = new Plot2D(width, height-50);
		testplot.showAxes(true);
		testplot.setLabels("dist. from next boat", "force on test boat");
		testplot.setTitle("D="+String.format("%.2f", D)+" k="+String.format("%.4f", k));
		testplot.boxed(true);
		testplot.lines(true);
		
		int testind = boats.size();
		double testdx = 0.1;
		int testpoints = 100;
		double[] testforce = new double[testpoints];
		double[] spacing = new double[testpoints];
		for(int i=0; i< testpoints; i++){
			spacing[i] = testdx*(i+1);
		}
		testplot.addX(spacing);
		
		int seriespoints = 12;
		for(int j=0; j<seriespoints; j++){
			//set up the cluster for these parameters
			k = k+.2;
			diag = 2*(D/dx/dx)+k;
			upper = -D/dx/dx + v/dx;
			lower = -D/dx/dx - v/dx;
			
			for(int h=0; h<testind; h++){
				removeBoat();
			}
			for(int h=0; h<testind; h++){
				addBoat();
			}
			shakeButton.doClick();
			shakeButton.doClick();
			shakeButton.doClick();
			calc();
			calc();
			int iter = 0;
			while(force[force.length-1] < force[force.length-2]*(0.95) && iter < 5){
				calc();
				iter++;
			}
			
			addBoat();
			
			//move the test boat and check force
			boats.get(testind).x = boats.get(testind-1).x - L - rc;
			for(int i=0; i< testpoints; i++){
				boats.get(testind).x -= testdx;
				
				setaF();
				solvec();
				solveForce();
				
				testforce[i] = force[testind] - force[testind-1];
			}
			
			//do something with the data
			if(j%4 == 0){
				testplot.addY(testforce, testplot.CIRCLE, Color.BLUE);
			}
			if(j%4==1){
				testplot.addY(testforce, testplot.CIRCLE, Color.GREEN);
			}
			if(j%4==2){
				testplot.addY(testforce, testplot.CIRCLE, Color.RED);
			}
			if(j%4==3){
				testplot.addY(testforce, testplot.CIRCLE, Color.BLACK);
			}
			//testplot.addY(testforce);
			
			removeBoat();
		}
		testplot.setXLimits(0.0, spacing[testpoints-1]+testdx);
		testplot.setYLimits(-0.025, 0.025);
		testplot.paint();
		PlotWindow pwind = new PlotWindow(testplot);
	}
	
	public void testNSeries(){
		Plot2D testplot = new Plot2D(width, height-50);
		testplot.showAxes(true);
		testplot.setLabels("dist. from next boat", "force on test boat");
		testplot.setTitle("D="+String.format("%.2f", D)+" k="+String.format("%.4f", k));
		testplot.boxed(true);
		testplot.lines(true);
		
		
		double testdx = 0.1;
		int testpoints = 100;
		double[] testforce = new double[testpoints];
		double[] spacing = new double[testpoints];
		double minf = 10000;
		double maxf = -minf;
		for(int i=0; i< testpoints; i++){
			spacing[i] = testdx*(i+1);
		}
		testplot.addX(spacing);
		
		int seriespoints = 10;
		for(int j=0; j<seriespoints; j++){
			
			addBoat();
			
			shakeButton.doClick();
			calc();
			calc();
			
			int testind = boats.size();
			addBoat();
			
			//move the test boat and check force
			boats.get(testind).x = boats.get(testind-1).x - L - rc;
			for(int i=0; i< testpoints; i++){
				boats.get(testind).x -= testdx;
				
				setaF();
				solvec();
				solveForce();
				
				testforce[i] = force[testind] - force[testind-1];
				minf = Math.min(minf, testforce[i]);
				maxf = Math.max(maxf, testforce[i]);
			}
			
			//do something with the data
			if(j%4 == 0){
				testplot.addY(testforce, testplot.CIRCLE, Color.BLUE);
			}
			if(j%4==1){
				testplot.addY(testforce, testplot.CIRCLE, Color.GREEN);
			}
			if(j%4==2){
				testplot.addY(testforce, testplot.CIRCLE, Color.RED);
			}
			if(j%4==3){
				testplot.addY(testforce, testplot.CIRCLE, Color.BLACK);
			}
			//testplot.addY(testforce);
			
			removeBoat();
		}
		testplot.setXLimits(0.0, spacing[testpoints-1]+testdx);
		testplot.setYLimits(minf-Math.min(Math.abs(minf), Math.abs(maxf))*0.05, maxf+Math.min(Math.abs(minf), Math.abs(maxf))*0.05);
		testplot.paint();
		PlotWindow pwind = new PlotWindow(testplot);
	}
	
	public void toInfinity(){
		
	}
	
	////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////

	public void addBoat(){
		Boat bn = new Boat(boats.size());
		if(bn.index == 0){
			bn.x = R-2;
			bn.v = v;
			force = new double[1];
			force[0] = 0.0;
		}else{
			bn.x = boats.get(boats.size()-1).x - (L+rc) - 0.1;
			bn.v = v;
			double[] tmp = force;
			force = new double[boats.size()+1];
			System.arraycopy(tmp, 0, force, 0, tmp.length);
			force[force.length-1] = force[force.length-2];
		}
		boats.add(bn);
		parPan.setN(boats.size());
		parPan.repaint();
	}

	public void removeBoat(){
		if(boats.isEmpty()){return;}
		boats.remove(boats.size()-1);
		force = new double[boats.size()];
		parPan.setN(boats.size());
		parPan.repaint();
	}

	public void writeImage(){
		if(plot != null){
			File imageFile = new File("clusterImages");
			if(!imageFile.exists()){
				imageFile.mkdir();
			}
			imageFile = new File("clusterImages/cluster"+String.valueOf(System.currentTimeMillis())+".png");	//for everything else
			try{
				ImageIO.write(plot, "png", imageFile);
				JOptionPane.showMessageDialog(frame, "successfully recorded image\nFile: "+imageFile.getName());
			}catch(Exception e){
				System.err.println("Error: " + e.getMessage());
				JOptionPane.showMessageDialog(frame, "file writing error");
			}
		}
	}

	public void setR(double newR){
		R = newR;
		nx = (int)(R*120);
		dx = R/(nx-1);

		x = new double[nx];
		c = new double[nx];
		y = new double[nx];
		q = new double[nx];
		aF = new double[nx];
		mc = new double[nx];
		md = new double[nx];
		for(int i=0; i<nx; i++){
			x[i] = i*dx;
			c[i] = 0.0;
			y[i] = 0.0;
			q[i] = 0.0;
			aF[i] = 0.0;
			mc[i] = 0.0;
			md[i] = 0.0;
		}
	}

	//now for some heavy math
	//this solves the almost tridiagonal system
	public void solvec(){
		if(boats.isEmpty()){
			for (int i = 0; i < nx; i++){
				c[i] = 0.0;
			}
			return;
		}

		double id = 0.0;
		//solve for y
		mc[0] = upper/(diag-1);
		md[0] = aF[0]/diag;
		for (int i = 1; i < nx; i++){
			id = 1/(diag - mc[i-1]*lower);
			if(i==nx-1){id = 1/(diag - upper*lower - mc[i-1]*lower);}
			mc[i] = upper*id;
			md[i] = (aF[i] - md[i-1]*lower)*id;
        }
        y[nx-1] = md[nx-1];
        for (int i=nx-2; i>=0; i--){
			y[i]=md[i]-mc[i]*y[i+1];
		}

		//solve for q
		mc[0] = upper/(diag-1);
		md[0] = 1.0/diag;
		for (int i = 1; i < nx; i++){
			id = 1/(diag - mc[i-1]*lower);
			if(i==nx-1){id = 1/(diag - upper*lower - mc[i-1]*lower);}
			mc[i] = upper*id;
			md[i] = (0.0 - md[i-1]*lower)*id;
        }
		md[nx-1] = (upper - md[nx-2]*lower)*id;
        q[nx-1] = md[nx-1];
        for (int i=nx-2; i>=0; i--){
			q[i]=md[i]-mc[i]*q[i+1];
		}

		//find the products
		double vy = y[0]+y[nx-1]*lower;
		double vq = q[0]+q[nx-1]*lower;

		//finally solve for c
		for (int i = 0; i < nx; i++){
			c[i] = y[i] - (vy/(1+vq))*q[i];
        }
	}
	
	//this solves the tridiagonal system 
	public void solvecInf(){
		
		
		
	}

	public void setaF(){
		for (int i = 0; i < nx; i++){
			aF[i] = 0.0;
        }
		for(Boat bn: boats){
			int idb = (int)Math.round((bn.x-rc)/dx);
			int idf = (int)Math.round((bn.x+rc)/dx);
			if(idb >= nx){idb = idb - nx;}
			if(idf >= nx){idf = idf - nx;}
			if(idb < 0){idb = idb + nx;}
			if(idf < 0){idf = idf + nx;}

			if(idf>idb){
				for(int i=idb; i<= idf; i++){
					aF[i] = a;
				}
			} else{
				for(int i=idb; i<nx; i++){
					aF[i] = a;
				}
				for(int i=0; i<= idf; i++){
					aF[i] = a;
				}
			}
		}
	}

	public void solveForce(){
		//first find the c difference for each boat
		for(Boat bn: boats){
			int idb = (int)Math.round((bn.x)/dx);
			int idf = (int)Math.round((bn.x+L)/dx);
			if(idb >= nx){idb = idb - nx;}
			if(idf >= nx){idf = idf - nx;}
			if(idb < 0){idb = idb + nx;}
			if(idf < 0){idf = idf + nx;}

			//the c difference
			//force[bn.index] = c[idf] - c[idb];
			//the tension difference
			//force[bn.index] = 22*(1/(b*b*c[idf]*c[idf]+1) - 1/(b*b*c[idb]*c[idb]+1));
			//the net force
			force[bn.index] = -v*vis + l*22*(1/(b*b*c[idf]*c[idf]+1) - 1/(b*b*c[idb]*c[idb]+1));
		}
	}

	/**
     *  This is a panel for displaying and modifying parameters
     */
	class PPanel extends JPanel {
		// there are lots of labels and fields
        JLabel NLabel, RLabel, mLabel, visLabel, LLabel, lLabel, rcLabel, DLabel, kLabel, aLabel, bLabel;
        public JTextField NField, RField, mField, visField, LField, lField, rcField, DField, kField, aField, bField;

		public PPanel(double[] params) {
			// there are lots of labels and fields
			NLabel = new JLabel("N = ");
            NLabel.setHorizontalAlignment(JLabel.TRAILING);
            RLabel = new JLabel("route length(cm) = ");
            RLabel.setHorizontalAlignment(JLabel.TRAILING);
            mLabel = new JLabel("mass(g) = ");
            mLabel.setHorizontalAlignment(JLabel.TRAILING);
            visLabel = new JLabel("viscosity constant(g/s) = ");
            visLabel.setHorizontalAlignment(JLabel.TRAILING);
            LLabel = new JLabel("length of boat(cm) = ");
            LLabel.setHorizontalAlignment(JLabel.TRAILING);
            lLabel = new JLabel("width of boat(cm) = ");
            lLabel.setHorizontalAlignment(JLabel.TRAILING);
            rcLabel = new JLabel("radius of camphor supply(cm) = ");
            rcLabel.setHorizontalAlignment(JLabel.TRAILING);
            DLabel = new JLabel("diffusion constant(cm*cm/s) = ");
            DLabel.setHorizontalAlignment(JLabel.TRAILING);
            kLabel = new JLabel("evaporation constant(1/s) = ");
            kLabel.setHorizontalAlignment(JLabel.TRAILING);
            aLabel = new JLabel("supply rate = ");
            aLabel.setHorizontalAlignment(JLabel.TRAILING);
            bLabel = new JLabel("beta = ");
            bLabel.setHorizontalAlignment(JLabel.TRAILING);

			// insert the initial parameters
			NField = new JTextField(String.valueOf(boats.size()), 6);
            NField.setEditable(true);
            RField = new JTextField(String.valueOf(params[0]), 6);
            RField.setEditable(true);
            mField = new JTextField(String.valueOf(params[1]), 6);
            mField.setEditable(true);
            visField = new JTextField(String.valueOf(params[2]), 6);
            visField.setEditable(true);
            LField = new JTextField(String.valueOf(params[3]), 6);
            LField.setEditable(true);
            lField = new JTextField(String.valueOf(params[4]), 6);
            lField.setEditable(true);
            rcField = new JTextField(String.valueOf(params[5]), 6);
            rcField.setEditable(true);
            DField = new JTextField(String.valueOf(params[6]), 6);
            DField.setEditable(true);
            kField = new JTextField(String.valueOf(params[7]), 6);
            kField.setEditable(true);
            aField = new JTextField(String.valueOf(params[8]), 6);
            aField.setEditable(true);
            bField = new JTextField(String.valueOf(params[9]), 6);
            bField.setEditable(true);

            // now set up the panel with the above items
            this.setLayout(new GridLayout(0,2));
			this.add(NLabel);
			this.add(NField);
            this.add(RLabel);
            this.add(RField);
            this.add(mLabel);
            this.add(mField);
            this.add(visLabel);
            this.add(visField);
            this.add(LLabel);
            this.add(LField);
            this.add(lLabel);
            this.add(lField);
            this.add(rcLabel);
            this.add(rcField);
            this.add(DLabel);
            this.add(DField);
            this.add(kLabel);
            this.add(kField);
            this.add(aLabel);
            this.add(aField);
            this.add(bLabel);
            this.add(bField);
		}

		public void newParameters(){
			try{
				setR(Double.valueOf(RField.getText()));
				m = Double.valueOf(mField.getText());
				vis = Double.valueOf(visField.getText());
				L = Double.valueOf(LField.getText());
				l = Double.valueOf(lField.getText());
				rc = Double.valueOf(rcField.getText());
				D = Double.valueOf(DField.getText());
				k = Double.valueOf(kField.getText());
				a = Double.valueOf(aField.getText());
				b = Double.valueOf(bField.getText());
			} catch(Exception expt){
				JOptionPane.showMessageDialog(this, "Couldn't read parameters. Try again.");
			}
		}

		public void setN(int n){
			NField.setText(String.valueOf(boats.size()));
		}
	}
	
	
}
