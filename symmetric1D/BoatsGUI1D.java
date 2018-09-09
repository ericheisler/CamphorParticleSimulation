import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import javax.swing.*;
import java.util.*;
import java.io.*;
import javax.imageio.ImageIO;

/**
 * This contains all of the model information and GUI components
 * @author Heisler
 */
class BoatsGUI1D implements ActionListener{
    // These are the parts of the GUI
	BoatsModel1D model, lyaponovmodel;
	int N;
	ArrayList<Boat1D> b;
	double[] c;
	double phasev, jamt1, jamt2, dt1, dt2;
	double[] density, dT, dc, force, flow;
	double[] velocity, frontspacing, backspacing, minspacing, avespacing, timesteps;
	double[][] allvelocity, allspacing, allposition, alllyapunov;

	BoatsWorker worker;

	JFrame frame;
	BPanel boatsPanel;
	PPanel paramPanel;
	JPanel cntrlPanel, rightPanel, bottomPanel;
	BufferedImage spaceTime, experiment, graphs, labelSpace, recordst;
	Plot2D p2d;
	int imageType;
	String[] imageNames;
	int simWidth, simHeight;
	double guiTimeStart, guiTimeEnd, modelTime;

	JMenuBar menuBar;
	JMenu fileMenu, modelMenu, windowMenu, simMenu, addOnMenu;
    JMenuItem infoItem, fRateItem, slowerItem, fastItem, sim1Item, sim2Item, sim3Item, sim4Item, sim5Item, sim6Item, equilItem, relaxItem;
	ButtonGroup windowGroup;
	JRadioButtonMenuItem bigWindowItem, smallWindowItem, recordImageItem;

    JButton runButton, timeButton, adButton, rmButton, homoButton, randButton;
	JButton paramButton, resetButton, stopButton, displayButton, backgButton, writeButton;

	JTextField backgField;

	boolean showc;
	long compytime;
	int frameRate;
	boolean computeFast, recordNext, recordData, computelyaponov;

	File sim1File, sim2File, imageFile;
	FileWriter writer, nondimWriter, spectWriter, corrWriter;

	/**
	 * The constructor initializes all the model parts and prepares
	 * everything to run.
	 */
	public BoatsGUI1D(){
		frameRate = 30;
		showc = false;
		computeFast = false;
		recordNext = false;
		recordData = false;
		computelyaponov = false;
		model = new BoatsModel1D();
		N = 0;
		b = model.getBoats();
		c = model.getc();
		phasev = 0.0;
		dt1 = 0.0;
		jamt1 = 0.0;

		worker = new BoatsWorker();

		simWidth = 800;
		simHeight = 600;
		imageType = 0;
		imageNames = new String[]{"space time", "experiment", "graphs", "c", "label space"};
		
		spaceTime = new BufferedImage(simWidth, simHeight, BufferedImage.TYPE_INT_ARGB);
		spaceTime.createGraphics().setColor(Color.WHITE);
		spaceTime.createGraphics().fillRect(0, 0, simWidth, simHeight);
		
		experiment = new BufferedImage(simWidth, simHeight, BufferedImage.TYPE_INT_ARGB);
		experiment.createGraphics().setColor(Color.WHITE);
		experiment.createGraphics().fillRect(0, 0, simWidth, simHeight);
		
		labelSpace = new BufferedImage(simWidth, simHeight, BufferedImage.TYPE_INT_ARGB);
		labelSpace.createGraphics().setColor(Color.WHITE);
		labelSpace.createGraphics().fillRect(0, 0, simWidth, simHeight);
		
		graphs = new BufferedImage(simWidth, simHeight, BufferedImage.TYPE_INT_ARGB);
		graphs.createGraphics().setColor(Color.WHITE);
		graphs.createGraphics().fillRect(0, 0, simWidth, simHeight);

		p2d= new Plot2D(simWidth, simHeight);
		p2d.setTitle("c(x)");
		p2d.setLabels("x", null);
		p2d.showAxes(true);
		p2d.lines(true);

		guiTimeStart = 0.0;
		guiTimeEnd = model.getInterval();
		modelTime = model.getTime();

		frame = new JFrame("boats simulation");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		//Create the menu bar.
        menuBar = new JMenuBar();
        menuBar.setOpaque(true);
		//create some menus
		fileMenu = new JMenu("file");
		modelMenu = new JMenu("model");
		windowMenu = new JMenu("window");
		simMenu = new JMenu("simulations");
		addOnMenu = new JMenu("add-ons");
		//create some menu items
        infoItem = new JMenuItem("model info");
		infoItem.addActionListener(this);

		sim1Item = new JMenuItem("fix N, vary vis(in time)");
		sim1Item.addActionListener(this);
		sim2Item = new JMenuItem("fix N, vary vis(in time), stat.");
		sim2Item.addActionListener(this);
		sim3Item = new JMenuItem("fix N vary vis and save images");
		sim3Item.addActionListener(this);
		sim4Item = new JMenuItem("fix N and vis, vary time");
		sim4Item.addActionListener(this);
		sim5Item = new JMenuItem("fix N and vis, lyap vs t");
		sim5Item.addActionListener(this);
		sim6Item = new JMenuItem("compute all the things!");
		sim6Item.addActionListener(this);
		
		equilItem = new JMenuItem("equilibrium");
		equilItem.addActionListener(this);
		relaxItem = new JMenuItem("constant v relaxor");
		relaxItem.addActionListener(this);

		fRateItem = new JMenuItem("frame rate");
		fRateItem.addActionListener(this);
		slowerItem = new JMenuItem("slow it down");
		slowerItem.addActionListener(this);
		fastItem = new JMenuItem("compute fast");
		fastItem.addActionListener(this);
		recordImageItem = new JRadioButtonMenuItem("start recording image");
		recordImageItem.setSelected(false);
		recordImageItem.addActionListener(this);
		windowGroup = new ButtonGroup();
		bigWindowItem = new JRadioButtonMenuItem("large window");
		bigWindowItem.setSelected(true);
		bigWindowItem.addActionListener(this);
		windowGroup.add(bigWindowItem);
		smallWindowItem = new JRadioButtonMenuItem("small window");
		smallWindowItem.setSelected(false);
		smallWindowItem.addActionListener(this);
		windowGroup.add(smallWindowItem);

		//put the menus together
		modelMenu.add(infoItem);
		windowMenu.add(fRateItem);
		windowMenu.add(slowerItem);
		windowMenu.add(fastItem);
		windowMenu.addSeparator();
		windowMenu.add(recordImageItem);
		windowMenu.addSeparator();
		windowMenu.add(bigWindowItem);
		windowMenu.add(smallWindowItem);
		simMenu.add(sim1Item);
		simMenu.add(sim2Item);
		simMenu.add(sim3Item);
		simMenu.add(sim4Item);
		simMenu.add(sim5Item);
		simMenu.add(sim6Item);
		addOnMenu.add(equilItem);
		addOnMenu.add(relaxItem);

		menuBar.add(fileMenu);
		menuBar.add(modelMenu);
		menuBar.add(windowMenu);
		menuBar.add(simMenu);
		menuBar.add(addOnMenu);

		//Set the menu bar
        frame.setJMenuBar(menuBar);

		//create some buttons
		runButton = new JButton("run");
		runButton.addActionListener(this);
		timeButton = new JButton("set T interval");
		timeButton.addActionListener(this);
		adButton = new JButton("add boat");
		adButton.addActionListener(this);
		rmButton = new JButton("remove boat");
		rmButton.addActionListener(this);
		homoButton = new JButton("homogenize");
		homoButton.addActionListener(this);
		randButton = new JButton("randomize");
		randButton.addActionListener(this);
		paramButton = new JButton("set parameters");
		paramButton.addActionListener(this);
		resetButton = new JButton("reset");
		resetButton.addActionListener(this);
		stopButton = new JButton("stop");
		stopButton.addActionListener(this);
		displayButton = new JButton("space time");
		displayButton.addActionListener(this);
		backgButton = new JButton("add constant c");
		backgButton.addActionListener(this);
		writeButton = new JButton("write to file");
		writeButton.addActionListener(this);
		
		backgField = new JTextField(String.valueOf(0.0), 6);
		backgField.setEditable(true);

		//make the boats panel
		boatsPanel = new BPanel();
		//make a control panel
		rightPanel = new JPanel(); //this holds the following panels
		cntrlPanel = new JPanel();
		paramPanel = new PPanel();
		bottomPanel = new JPanel();

		cntrlPanel.setLayout(new GridLayout(0,2));
		cntrlPanel.add(runButton);
		cntrlPanel.add(timeButton);
		cntrlPanel.add(adButton);
        cntrlPanel.add(rmButton);
		cntrlPanel.add(homoButton);
		cntrlPanel.add(randButton);
		cntrlPanel.add(displayButton);
		cntrlPanel.add(resetButton);
		cntrlPanel.add(writeButton);
		cntrlPanel.add(stopButton);
		cntrlPanel.add(backgButton);
		cntrlPanel.add(paramButton);
		cntrlPanel.add(backgField);

		rightPanel.setLayout(new BorderLayout());
		rightPanel.add(cntrlPanel, BorderLayout.NORTH);
		rightPanel.add(paramPanel, BorderLayout.CENTER);
		rightPanel.add(bottomPanel, BorderLayout.SOUTH);

		// for the panels
		boatsPanel.setOpaque(true);
		cntrlPanel.setOpaque(true);
		rightPanel.setOpaque(true);
		bottomPanel.setOpaque(true);

		Container pane = frame.getContentPane();
		pane.setLayout(new BorderLayout());
		pane.add(boatsPanel, BorderLayout.CENTER);
		pane.add(rightPanel, BorderLayout.LINE_END);

		pane.validate();
		//Display the window.
        frame.pack();
        frame.setVisible(true);
	}

	// for the actionlistener
	public void actionPerformed(ActionEvent e){
		if (e.getSource() == runButton){
			if(N > 0 && !worker.getState().equals(SwingWorker.StateValue.STARTED)){
				worker = new BoatsWorker();
				worker.execute();
			}
		} else if (e.getSource() == timeButton){
			double t = Double.valueOf((String)JOptionPane.showInputDialog(frame, "Input new time interval"));
			model.setInterval(t);
			guiTimeStart = model.getTime();
			guiTimeEnd = guiTimeStart+t;
			boatsPanel.repaint();
		} else if (e.getSource() == adButton){
			if(worker != null){
				worker.cancel(true);
			}
			N = model.addBoat();
			b = model.getBoats();
			paramPanel.NField.setText(String.valueOf(N));
			paramPanel.repaint();
			boatsPanel.repaint();
		} else if (e.getSource() == rmButton){
			if(worker != null){
				worker.cancel(true);
			}
			if(N>0){
				N = model.removeBoat();
				b = model.getBoats();
			}
			paramPanel.NField.setText(String.valueOf(N));
			paramPanel.repaint();
			boatsPanel.repaint();
		} else if (e.getSource() == homoButton){
			model.homogenize();
			boatsPanel.repaint();
		} else if (e.getSource() == randButton){
			model.randomize();
			boatsPanel.repaint();
		} else if (e.getSource() == paramButton){
			paramPanel.newParameters();
		} else if (e.getSource() == resetButton){
			reset();
			paramPanel.NField.setText(String.valueOf(N));
			paramPanel.repaint();
		} else if (e.getSource() == stopButton){
			if(worker != null){
				worker.cancel(true);
			}
		} else if (e.getSource() == writeButton){
			if(worker != null){
				worker.cancel(true);
			}
			computeData();
			writeData();
		} else if (e.getSource() == displayButton){
			imageType = (imageType+1)%imageNames.length;
			displayButton.setText(imageNames[imageType]);
			if(imageType == 4 || imageType == 5){
				recordData = true;
			}else{
				recordData = false;
			}
			frame.repaint();
		} else if (e.getSource() == backgButton){
			model.kickABoat(0,1000.0);
			try{
				model.changeBackground(Double.valueOf(backgField.getText()));
			} catch(Exception expt){
				JOptionPane.showMessageDialog(paramPanel, "Couldn't read parameter. Try again.");
			}
		} else if (e.getSource() == infoItem){
			JOptionPane.showMessageDialog(frame, "Hold on a moment, I will explain.");
		} else if (e.getSource() == fRateItem){
			frameRate = Integer.valueOf((String)JOptionPane.showInputDialog(frame, "Input an integer 10-60\n (this is the frame rate in model time, NOT real time)"));
		} else if (e.getSource() == fastItem){
			computeFast = !computeFast;
		} else if (e.getSource() == recordImageItem){
			if(!worker.getState().equals(SwingWorker.StateValue.STARTED)){
				recordNext = !recordNext;
				if(recordNext){
					recordImageItem.setText("stop recording image");
					
				}else{
					recordImageItem.setText("start recording image");
					double[] pars = model. getParameters();
					writeImage(recordst, "images/N"+String.valueOf(b.size())+"vis"+String.valueOf((int)(10000*pars[2]))+"K"+String.valueOf((int)pars[7])+"T"+String.valueOf((int)model.getTime())+".png");
					recordst = null;
				}
			}
			
		} else if (e.getSource() == bigWindowItem){
			boatsPanel.setSize(900, 800);
			simHeight = 600;
			reset();
			N = 0;
			guiTimeStart = model.getTime();
			guiTimeEnd = guiTimeStart+model.getInterval();
			frame.pack();
			frame.repaint();
		} else if (e.getSource() == smallWindowItem){
			boatsPanel.setSize(900, 600);
			simHeight = 400;
			reset();
			N = 0;
			guiTimeStart = model.getTime();
			guiTimeEnd = guiTimeStart+model.getInterval();
			frame.pack();
			frame.repaint();
		} else if (e.getSource() == sim1Item){
			runSim1();
		} else if (e.getSource() == sim2Item){
			runSim2();
		} else if (e.getSource() == sim3Item){
			runSim3();
		} else if (e.getSource() == sim4Item){
			runSim4();
		} else if (e.getSource() == sim5Item){
			runSim5();
		} else if (e.getSource() == sim6Item){
			runSim6();
		} else if (e.getSource() == equilItem){
			//Equilibrator equil = new Equilibrator(model.getParameters());
		} else if (e.getSource() == relaxItem){
			//Relaxor laxor = new Relaxor(model.getParameters());
		} 
	}
	
	/////////////////////////////////////////////////////////////////////////////////////////
	// beginning of simulations //
	/////////////////////////////////////////////////////////////////////////////////////////

	/**
	 * 
	 */
	public void runSim1(){
				
	}
	
	
	/**
	 * 
	 */
	public void runSim2(){
		
	}
	
	/**
	 * 
	 */
	public void runSim3(){
		
	}
	
	/**
	 * 
	 */
	public void runSim4(){
		
	}
	
	/**
	 * 
	 */
	public void runSim5(){
		
	}
	
	/**
	 * 
	 */
	public void runSim6(){
		
	}
	
	//////////////////////////////////////////////////////////////////////////////////////
	// end of simulations //
	//////////////////////////////////////////////////////////////////////////////////////

	public void computeData(){
		double[] dims = model.getParameters();
		double R = dims[0];
		double l = dims[3];
		double L = dims[4];
		double rc = dims[5];
		double beta = dims[9];
		double vis = dims[2];
		dims = null;
		double dx = model.getDisc()[4];
		int nx = (int)model.getDisc()[5];

		int indb, indf;
		double Tb, Tf;
		density = new double[b.size()];
		dT = new double[b.size()];
		dc = new double[b.size()];
		force = new double[b.size()];
		flow = new double[b.size()];

		int adjdiam = (int)Math.round(L/dx);
		adjdiam += adjdiam%2;
		for(Boat1D bn: b){
			indb = (int)Math.round(bn.x/dx)-adjdiam/2;		// find the index for the back
			indf = indb+adjdiam;		// find the index for the front
			if(indb >= nx){ indb -= nx; }	// for periodic boundary
			if(indf >= nx){ indf -= nx; }	//
			if(indb < 0){ indb += nx; }	// for periodic boundary
			if(indf < 0){ indf += nx; }	//
			Tb = 22/((beta*c[indb])*(beta*c[indb]) + 1); //surface tension at front
			Tf = 22/((beta*c[indf])*(beta*c[indf]) + 1); //surface tension at back

			density[bn.index] = 2/(bn.back+bn.front+2*(L+rc));
			dT[bn.index] = Tf-Tb;
			dc[bn.index] = c[indb] - c[indf];
			force[bn.index] = (Tf-Tb)*l - vis*bn.v;
			flow[bn.index] = 60*bn.v*density[bn.index];
		}
	}
	
	// writes data to file
	public void writeData(){
		// nothing here
	}
	
	public void writeImage(BufferedImage img, String name){
		if(img != null){
			imageFile = new File("images");
			if(!imageFile.exists()){
				imageFile.mkdir();
			}
			imageFile = new File(name);	//for everything else
			try{
				ImageIO.write(img, "bmp", imageFile);
				//JOptionPane.showMessageDialog(frame, "successfully recorded image\nFile: "+imageFile.getName());
			}catch(Exception e){
				System.err.println("Error: " + e.getMessage());
				JOptionPane.showMessageDialog(frame, "file writing error");
			}
		}
	}
	
	public void reset(){
		worker.cancel(true);
		worker = new BoatsWorker();
		model.reset();
		b = model.getBoats();
		c = model.getc();
		N = 0;
		guiTimeStart = model.getTime();
		guiTimeEnd = guiTimeStart+model.getInterval();
		modelTime = model.getTime();
		phasev = 0.0;
		spaceTime = null;
		spaceTime = new BufferedImage(simWidth, simHeight, BufferedImage.TYPE_INT_ARGB);
		spaceTime.createGraphics().fillRect(0, 0, simWidth, simHeight);
		recordst = null;
		experiment = null;
		experiment = new BufferedImage(simWidth, simHeight, BufferedImage.TYPE_INT_ARGB);
		experiment.createGraphics().fillRect(0, 0, simWidth, simHeight);
		graphs = null;
		graphs = new BufferedImage(simWidth, simHeight, BufferedImage.TYPE_INT_ARGB);
		graphs.createGraphics().fillRect(0, 0, simWidth, simHeight);
		labelSpace = null;
		labelSpace = new BufferedImage(simWidth, simHeight, BufferedImage.TYPE_INT_ARGB);
		labelSpace.createGraphics().fillRect(0, 0, simWidth, simHeight);
		p2d= new Plot2D(simWidth, simHeight);
		p2d.setTitle("c(x)");
		p2d.setLabels("x", null);
		p2d.setXLimits(0.0, 45.5);
		p2d.showAxes(true);
		p2d.lines(true);

		boatsPanel.repaint();
	}

	//////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////
	/**
	 * This is the worker thread which runs the model calculations
	 * in the background and publishes the results
	 */
	private class BoatsWorker extends SwingWorker<String, Object[]> {
		boolean drawn;
		BufferedImage recordtemp;
		
		@Override
        protected String doInBackground() {
            double[] disc = model.getDisc();
			int totalFrames = 1;
			int stepsPerChunk = (int)disc[2];
			int extraSteps = 0;
			if(frameRate > 0){
				totalFrames = (int)(disc[0]*frameRate);
				stepsPerChunk = (int)(disc[2] / totalFrames);
				extraSteps = ((int)disc[2])%(stepsPerChunk*totalFrames);
				while(extraSteps >= stepsPerChunk){
					extraSteps -= stepsPerChunk;
					totalFrames++;
				}
			}
			compytime = System.currentTimeMillis();
			double modeltime = model.getTime();
			double stTime = 0.0;
			guiTimeStart = model.getTime();
			guiTimeEnd = guiTimeStart+model.getInterval();
			Color[] cColor = new Color[(int)disc[5]];
			
			ArrayList<Boat1D> bgb = new ArrayList<Boat1D>(N);
			double[] bgc = new double[(int)disc[5]];
			
			//for an extra fast computation without the fancy stuff
			if(computeFast){
				double delta = 1e-3;
				double tmpdc;
				double[] tmpx = new double[N];
				double[] tmpv = new double[N];
				int stepsPerStep = (int)(5.0*disc[2]/disc[0]);
				int lyapSteps = (int)(disc[2]/stepsPerStep);
				double[] deltax = new double[N];
				double[] deltav = new double[N];
				double[] deltac = new double[N];
				if(recordData){
					velocity = new double[(int)(disc[2]/4)];
					frontspacing = new double[velocity.length];
					backspacing = new double[velocity.length];
					minspacing = new double[velocity.length];
					avespacing = new double[velocity.length];
					timesteps = new double[velocity.length];
					allvelocity = new double[velocity.length][N];
					allspacing = new double[velocity.length][N];
					allposition = new double[velocity.length][N];
					if(computelyaponov){
						lyaponovmodel = model.makeACopy();
						alllyapunov = new double[lyapSteps][N];
						lyaponovmodel.copyc(model.getc());
						lyaponovmodel.copyBoats(model.getBoats());
						for(int i=0; i<N; i++){
							deltax[i] = delta/1.41421356;
							deltav[i] = delta/1.41421356;
							deltac[i] = delta;
							lyaponovmodel.getBoats().get(i).x += deltax[i];
							lyaponovmodel.getBoats().get(i).v += deltav[i];
						}
						for(int i=0; i<N; i++){
							lyaponovmodel.update(i);
						}
					}
				}
				int vsind = 0;
				int lsind = 0;
				double routeR = model.getDisc()[3];
				
				// do all the time steps and store data, no graphical stuff
				for(int i=0; i<(int)disc[2]; i++){
					model.timeStep();
					
					if(recordData){
						//record data every 4 steps
						if(i%4 == 0 && vsind < velocity.length){
							bgb = model.getBoats();
							
							velocity[vsind] = bgb.get(0).v;
							frontspacing[vsind] = bgb.get(0).front;
							backspacing[vsind] = bgb.get(0).back;
							minspacing[vsind] = Math.min(frontspacing[vsind], backspacing[vsind]);
							avespacing[vsind] = (frontspacing[vsind] + backspacing[vsind])/2;
							for(int j=0; j<N; j++){
								allvelocity[vsind][j] = bgb.get(j).v;
								allspacing[vsind][j] = bgb.get(j).front;
								allposition[vsind][j] = bgb.get(j).x;
							}
							
							timesteps[vsind] = i*disc[1];
							vsind++;
						}
						//compute lyaponov stuff
						if(computelyaponov && i%stepsPerStep == 0 && lsind < alllyapunov.length){
							if(lsind > 0){
								for(int j=0; j<N; j++){
									tmpx[j] = lyaponovmodel.getBoats().get(j).front;
									tmpv[j] = lyaponovmodel.getBoats().get(j).v;
								}
								
								for(int j=0; j<stepsPerStep; j++){
									lyaponovmodel.timeStep();
								}
								
								for(int j=0; j<N; j++){
									deltax[j] = lyaponovmodel.getBoats().get(j).front - tmpx[j];
									//if(deltax[j] > routeR/2){ deltax[j] -= routeR; }
									//if(deltax[j] < -routeR/2){ deltax[j] += routeR; }
									deltav[j] = lyaponovmodel.getBoats().get(j).v - tmpv[j];
									tmpdc = Math.sqrt(deltax[j]*deltax[j] + deltav[j]*deltav[j] + 1e-99);
									alllyapunov[lsind][j] = Math.log(tmpdc/deltac[j]);
									deltac[j] = tmpdc;
								}
								
								lyaponovmodel.copyc(model.getc());
								lyaponovmodel.copyBoats(model.getBoats());
								
								for(int j=0; j<N; j++){
									lyaponovmodel.getBoats().get(j).x -= delta*deltax[j]/deltac[j];
									if(lyaponovmodel.getBoats().get(j).x > routeR){
										lyaponovmodel.getBoats().get(j).x -= routeR;
									}
									if(lyaponovmodel.getBoats().get(j).x < 0){
										lyaponovmodel.getBoats().get(j).x += routeR;
									}
									lyaponovmodel.getBoats().get(j).v += delta*deltav[j]/deltac[j];
								}
								for(int j=0; j<N; j++){
									lyaponovmodel.update(j);
								}
							}
							lsind++;
						}
					}
				}
				stTime += disc[2]*disc[1];
				bgb = model.getBoats();
				bgc = model.getc();
				Object[] vxcst = {bgb, bgc, new Double(model.getTime()), null};
				publish(vxcst);
				
				return "finished fast";
			}
			
			if(recordNext){
				int length = (int)Math.max(simHeight, simHeight*model.getInterval()/10);
				recordtemp = new BufferedImage(simWidth, length, BufferedImage.TYPE_INT_ARGB);
				Graphics2D rtg = recordtemp.createGraphics();
				rtg.setColor(Color.WHITE);
				rtg.fillRect(0, 0, simWidth, length);
			}
			drawn = false;
			
			float minc = 10;
			float maxc = 1;
			double maxfc = 5;
			
			BufferedImage bgst;
			int[] xpix;
			int tpix;
			bgst = new BufferedImage(simWidth, simHeight, BufferedImage.TYPE_INT_ARGB);
			bgst.createGraphics().setColor(Color.WHITE);
			bgst.createGraphics().fillRect(0, 0, simWidth, simHeight);
			bgst.createGraphics().setColor(Color.BLACK);
			
			xpix = new int[N];
			
			if(recordData){
				velocity = new double[(int)(disc[2]/4)];
				frontspacing = new double[velocity.length];
				backspacing = new double[velocity.length];
				minspacing = new double[velocity.length];
				avespacing = new double[velocity.length];
				timesteps = new double[velocity.length];
				allvelocity = new double[velocity.length][N];
				allspacing = new double[velocity.length][N];
				allposition = new double[velocity.length][N];
			}
			int vsind = 0;
			double delta = 0.0000001;
			double deltax, deltav, deltac;
			
			//start performing time steps until reaching T
			for(int j=0; j<totalFrames+1; j++){
				if(isCancelled()){
					System.out.println("got cancelled!");
					return "cancelled";
				}
				if(j == totalFrames){
					stepsPerChunk = extraSteps;
				}
				for(int i=0; i<stepsPerChunk; i++){
					model.timeStep();	//do several steps for each frame
					// but remember to record the data for each step
					bgb = model.getBoats();
					bgc = model.getc();
					tpix = (int)((bgst.getHeight()-1)*stTime/(guiTimeEnd-guiTimeStart)+1);
					
					//record an image for saving
					if(recordNext){
						int rtpix = (int)((recordtemp.getHeight()-1)*stTime/(guiTimeEnd-guiTimeStart)+1);
						xpix[0] = (int)((recordtemp.getWidth()-1)*bgb.get(0).x/disc[3])+1;
						recordtemp.setRGB(xpix[0], recordtemp.getHeight()-rtpix, Color.BLUE.getRGB());
						for(int k=1; k<N; k++){
							xpix[k] = (int)((recordtemp.getWidth()-1)*bgb.get(k).x/disc[3])+1;
							recordtemp.setRGB(xpix[k], recordtemp.getHeight()-rtpix, Color.RED.getRGB());
						}
					}
					
					//draw the boat positions
					xpix[0] = (int)((bgst.getWidth()-1)*bgb.get(0).x/disc[3])+1;
					bgst.setRGB(xpix[0], bgst.getHeight()-tpix, Color.BLUE.getRGB());
					for(int k=1; k<N; k++){
						xpix[k] = (int)((bgst.getWidth()-1)*bgb.get(k).x/disc[3])+1;
						bgst.setRGB(xpix[k], bgst.getHeight()-tpix, Color.RED.getRGB());
					}
					
					//store velocity and spacing data for the first boat
					if(recordData){
						if((j*stepsPerChunk + i)%4 == 0 && vsind < velocity.length){
							velocity[vsind] = bgb.get(0).v;
							frontspacing[vsind] = bgb.get(0).front;
							backspacing[vsind] = bgb.get(0).back;
							minspacing[vsind] = Math.min(frontspacing[vsind], backspacing[vsind]);
							avespacing[vsind] = frontspacing[vsind];
							for(int l=0; l<N; l++){
								allvelocity[vsind][l] = bgb.get(l).v;
								allspacing[vsind][l] = bgb.get(l).front;
								allposition[vsind][l] = bgb.get(l).x;
							}
							timesteps[vsind] = stTime;
							vsind++;
						}
					}
										
					stTime += disc[1];
				}
				bgb = model.getBoats();
				bgc = model.getc();
				Object[] vxcst = {bgb, bgc, new Double(model.getTime()), bgst};
				publish(vxcst);
				
			}
            return "finished"; //just because they said to do it
        }
		
        @Override
        protected void process(java.util.List<Object[]> vxcst) {
			Object[] billy = vxcst.get(vxcst.size()-1);
			b = (ArrayList<Boat1D>) billy[0];
			c = (double[]) billy[1];
			modelTime = ((Double) billy[2]).doubleValue();
			if(!computeFast){
				spaceTime.getRaster().setRect(((BufferedImage) billy[3]).getRaster());
				//spaceTime = (BufferedImage) billy[3];
			}
			
			//            b = (Vector<Boat>)((Vector<Boat>) billy[0]).clone();
			//			c = (Vector<Double>)((Vector<Double>) billy[1]).clone();
			//			modelTime = ((Double) billy[2]).doubleValue();
			//			spaceTime = ((BufferedImage) billy[3]).getSubimage(0,0,((BufferedImage)billy[3]).getWidth(),((BufferedImage)billy[3]).getHeight());
			
			double maxv = -100.0;
			double minv = 100.0;
			for(Boat1D bn: b){
				maxv = Math.max(maxv, bn.v);
				minv = Math.min(minv, bn.v);
			}
			boatsPanel.paintImmediately(0,0,boatsPanel.getWidth(),boatsPanel.getHeight());
        }
		
		@Override
		protected void done() {
			try{
				get();
			}catch(Exception ignore){}
			
			// add the recordtemp to the recordst image
			if(recordNext && !drawn){
				drawn = true;
				int length = (int)Math.max(simHeight, simHeight*model.getInterval()/10);
				if(recordst != null){
					// expand recordst
					BufferedImage bigger = new BufferedImage(simWidth, length+recordst.getHeight(), BufferedImage.TYPE_INT_ARGB);
					Graphics2D bigg = bigger.createGraphics();
					bigg.drawImage(recordst, 0, length, null);
					bigg.drawImage(recordtemp, 0, 0, null);
					recordst = bigger;
				}else if(recordtemp != null){
					recordst = new BufferedImage(simWidth, length, BufferedImage.TYPE_INT_ARGB);
					recordst.createGraphics().drawImage(recordtemp, 0, 0, null);
				}else{
					drawn = false;
				}
				
			}
			
			boatsPanel.paintImmediately(0,0,boatsPanel.getWidth(),boatsPanel.getHeight());
		}
		
	}
	
//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	/**
     * This is the graphical image of the boat experiment
     * NOTE: no math is done through this class, only drawing
     */
	class BPanel extends JPanel {

		int width;
		int height;

		public BPanel() {
			width = 900;
			height = 850;
		}

		public Dimension getPreferredSize() {
			return new Dimension(width, height);
		}

		public void setSize(int w, int h){
			width = w;
			height = h;
		}

		// ye olde painting method
		public void paintComponent(Graphics g) {
			super.paintComponent(g);
			//draw things on an image
			BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
			Graphics ig = image.getGraphics();

			// get some model dimensions
			double[] dims = model.getParameters();
			double R = dims[0];
			double L = dims[3];
			double rc = dims[4];

			int esimHeight = (int)(L/R*simWidth) + 2;
			int simposx = 20;
			int simposy = height-70;

			int picposx = simposx;
			int picposy = 20;

			//normalize c for drawing
			float minc = Float.MAX_VALUE-1;
			float maxc = -minc;
			for(int i=0; i<c.length; i++){
				if(c[i]>maxc)
					maxc = (float)c[i];
				if(c[i]<minc)
					minc = (float)c[i];
			}
			if(minc>=maxc){
				minc = 0;
			}
			//and the bulk
			float minvc = Float.MAX_VALUE-1;
			float maxvc = -minvc;
			for(int i=0; i<c.length; i++){
				if(model.volc[i]>maxvc)
					maxvc = (float)model.volc[i];
				if(model.volc[i]<minvc)
					minvc = (float)model.volc[i];
			}
			if(minvc>=maxvc){
				minvc = 0;
			}
			//normalize v for drawing
			float minv = Float.MAX_VALUE-1;
			float maxv = -minv;
			for(Boat1D bn : b){
				maxv = Math.max(maxv, (float)bn.v);
				minv = Math.min(minv, (float)bn.v);
			}
			if(minv>=maxv){
				minv = 0;
			}

			//draw a linear version of the experiment
			ig.setColor(Color.WHITE);
			ig.drawRect(simposx, simposy, simWidth, esimHeight);
			ig.setColor(Color.BLACK);
			ig.drawRect(simposx, simposy, simWidth, esimHeight);
			
			//draw the camphor concentration
			ig.setColor(Color.WHITE);
			ig.fillRect(simposx, simposy-110, simWidth, 100);
			ig.setColor(Color.BLACK);
			ig.drawRect(simposx, simposy-110, simWidth, 100);
			// the surface
			ig.setColor(Color.RED);
			int cpix = 0;
			for(int i=0; i<c.length; i++){
				cpix = simposy - 10 - (int)(100*((c[i]-minc)/(maxc-minc)));
				ig.drawLine((int)(simposx+i*800/c.length), cpix, (int)(simposx+i*800/c.length), cpix);
			}
			ig.drawString("camphor concentration, threshold = "+String.format("%.4f",1/dims[11]), simposx+100, simposy-115);
			ig.drawString(String.format("%.4f",maxc), simWidth+simposx+5, simposy-100);
			ig.drawString(String.format("%.4f",minc), simWidth+simposx+5, simposy-25);
			// and the bulk
			ig.setColor(new Color(0, 100, 0));
			int vcpix = 0;
			for(int i=0; i<model.nx; i++){
				vcpix = simposy - 10 - (int)(100*((model.volc[i]-minvc)/(maxvc-minvc)));
				ig.drawLine((int)(simposx+i*800/c.length), vcpix, (int)(simposx+i*800/c.length), vcpix);
			}
			ig.drawString(String.format("%.4f",maxvc), simWidth+simposx+5, simposy-85);
			ig.drawString(String.format("%.4f",minvc), simWidth+simposx+5, simposy-10);

			//draw the velocity
			ig.setColor(Color.BLUE);
			int vpix = 0;
			for(Boat1D bn: b){
				vpix = simposy - 10 - (int)(100*((bn.v-minv)/(maxv-minv)));
				ig.fillOval(simposx+(int)(bn.x/R * simWidth), vpix-2, 5, 5);
			}
			ig.drawString("velocity(cm/s)", simposx+450, simposy-115);
			ig.drawString(String.format("%.4f",maxv), simWidth+simposx+5, simposy-70);
			ig.drawString(String.format("%.4f",minv), simWidth+simposx+5, simposy-40);

			//draw the boats
			ig.setColor(Color.BLACK);
			int csize = (int)(2*rc/R*simWidth);
			int bsize = (int)(L/R*simWidth);
			for(Boat1D bn: b){
				int bpos = (int) ((bn.x-L/2)/R * simWidth);
				int cpos = (int) ((bn.x-rc)/R * simWidth);
				ig.setColor(Color.BLACK);
				ig.fillOval(bpos+simposx, simposy+(int)(esimHeight/2 - bsize/2), bsize, bsize);
				ig.setColor(Color.WHITE);
				ig.drawOval(cpos+simposx, simposy+(int)(esimHeight/2 - csize/2), csize, csize);
				if(bn.index == 0){
					ig.setColor(Color.RED);
					ig.fillOval(bpos+simposx, simposy+(int)(esimHeight/2 - bsize/2), bsize, bsize);
					ig.setColor(Color.BLACK);
				}
			}

			//draw the big picture
			if(imageType == 0){
				ig.drawImage(spaceTime, picposx, picposy, this);
			}else if(imageType == 1){
				drawExperiment();
				ig.drawImage(experiment, picposx, picposy, this);
			}else if(imageType == 2){
				drawGraphs();
				ig.drawImage(graphs, picposx, picposy, this);
			}else if(imageType == 3){
				drawC();
				ig.drawImage(p2d, picposx, picposy, this);
			}else if(imageType == 4){
				drawLabelSpace();
				ig.drawImage(labelSpace, picposx, picposy, this);
			}

			//draw some useful stuff
			ig.setColor(Color.BLACK);
			ig.drawString("x(cm)", (int)(simWidth/2+simposx), spaceTime.getHeight()+30);

			ig.drawString("time(s)", simWidth+simposx+5, (int)(spaceTime.getHeight()/2+picposy));
			ig.drawString("t = "+String.format("%.2f",modelTime), simWidth+simposx+5, (int)(spaceTime.getHeight()/2+picposy + 20));
			ig.drawString(String.format("%.2f",guiTimeEnd), simWidth+simposx+5, picposy+10);
			ig.drawString(String.format("%.2f",guiTimeStart), simWidth+simposx+5, spaceTime.getHeight()+picposy);

			ig.drawString("N = "+String.valueOf(N), simposx+simWidth+10, simposy+10);

			g.drawImage(image, 0, 0, null);
		}
		
		public void drawExperiment(){
			double[] dims = model.getParameters();
			double R = dims[0];
			double L = dims[3];
			double rc = dims[4];
			Graphics2D exg = experiment.createGraphics();
			int exDia = experiment.getHeight()-40;
			double scale = (exDia)/(R/Math.PI + L);
			int boatDia = (int)(L*scale);
			exg.setColor(Color.WHITE);
			exg.fillRect(0, 0, experiment.getWidth(), experiment.getHeight());
			exg.setColor(Color.BLACK);
			exg.drawOval(20, 20, exDia, exDia);
			exg.drawOval(20+2*boatDia, 20+2*boatDia, exDia-4*boatDia, exDia-4*boatDia);
			exg.drawLine(50+exDia, 20, 50+exDia, 20+exDia);
			exg.drawLine(40+exDia, 20, 60+exDia, 20);
			exg.drawLine(40+exDia, 20+exDia, 60+exDia, 20+exDia);
			exg.drawString(String.format("%.2f", R/Math.PI)+"cm", 55+exDia, 10+exDia/2);
			for(Boat1D bn: b){
				int bxpos = (int) (20+(Math.cos(bn.x/R*2*Math.PI)+1)*(exDia-2*boatDia)/2 +boatDia/2);
				int bypos = (int) (20+(Math.sin(bn.x/R*2*Math.PI)+1)*(exDia-2*boatDia)/2 +boatDia/2);
				exg.fillOval(bxpos, bypos , boatDia, boatDia);
			}
		}
		
		public void drawGraphs(){
			Graphics2D grg = graphs.createGraphics();
			int graphHeight = (int)(graphs.getHeight()/5 - 15-1);
			int graphWidth = 800;
			int posy = 15;
			double R = model.getParameters()[0];
			double min = 0.0;
			double max = 1.0;
			double[] xval = new double[b.size()];
			for(Boat1D bn: b){
				xval[bn.index] = bn.x;
			}

			computeData();
			
			//clear the image
			grg.setColor(Color.WHITE);
			grg.fillRect(0, 0, graphs.getWidth(), graphs.getHeight());


			//draw density (average front and back)
			//normalize density for drawing
			min = Double.MAX_VALUE-1;
			max = -min;
			for(Boat1D bn : b){
				max = Math.max(max, density[bn.index]);
				min = Math.min(min, density[bn.index]);
			}
			if(min>=max){
				min = 0.0;
			}
			grg.setColor(Color.BLACK);
			grg.drawRect(0, posy, graphWidth-1, graphHeight);
			grg.setColor(Color.RED);
			int dpix = 0;
			for(Boat1D bn: b){
				dpix = posy +graphHeight - (int)(graphHeight*((density[bn.index]-min)/(max-min)));
				grg.fillOval((int)(bn.x/R * graphWidth), dpix-2, 5, 5);
			}
			grg.drawString("density(boats/cm)", 50, posy - 2);
			grg.drawString("max = "+String.format("%.2f",max), 400, posy - 2);
			grg.drawString("min = "+String.format("%.2f",min), 600, posy - 2);
			
			//draw change in surface tension
			posy += graphHeight+15;
			
			//normalize dT for drawing
			min = Double.MAX_VALUE-1;
			max = -min;
			for(Boat1D bn : b){
				max = Math.max(max, dT[bn.index]);
				min = Math.min(min, dT[bn.index]);
			}
			if(min>=max){
				min = 0.0;
			}
			grg.setColor(Color.BLACK);
			grg.drawRect(0, posy, graphWidth-1, graphHeight);
			grg.setColor(Color.RED);
			int dTpix = 0;
			for(Boat1D bn: b){
				dTpix = posy +graphHeight - (int)(graphHeight*((dT[bn.index]-min)/(max-min)));
				grg.fillOval((int)(bn.x/R * graphWidth), dTpix-2, 5, 5);
			}
			if(min < 0 && max > 0){
				int zeroline = posy +graphHeight - (int)(graphHeight*((0-min)/(max-min)));
				grg.drawLine(0, zeroline, graphWidth, zeroline);
			}
			grg.drawString("surface tension difference(g/s/s) front-back", 50, posy - 2);
			grg.drawString("max = "+String.format("%.2f",max), 400, posy - 2);
			grg.drawString("min = "+String.format("%.2f",min), 600, posy - 2);
			
			//draw change in c
			posy += graphHeight+15;
			//normalize dc for drawing
			min = Double.MAX_VALUE-1;
			max = -min;
			for(Boat1D bn : b){
				max = Math.max(max, dc[bn.index]);
				min = Math.min(min, dc[bn.index]);
			}
			if(min>=max){
				min = 0.0;
			}
			grg.setColor(Color.BLACK);
			grg.drawRect(0, posy, graphWidth-1, graphHeight);
			grg.setColor(Color.RED);
			int dcpix = 0;
			for(Boat1D bn: b){
				dcpix = posy +graphHeight - (int)(graphHeight*((dc[bn.index]-min)/(max-min)));
				grg.fillOval((int)(bn.x/R * graphWidth), dcpix-2, 5, 5);
			}
			if(min < 0 && max > 0){
				int zeroline = posy +graphHeight - (int)(graphHeight*((0-min)/(max-min)));
				grg.drawLine(0, zeroline, graphWidth, zeroline);
			}
			grg.drawString("camphor difference back-front", 50, posy - 2);
			grg.drawString("max = "+String.format("%.2f",max), 400, posy - 2);
			grg.drawString("min = "+String.format("%.2f",min), 600, posy - 2);
			
			//draw net force
			posy += graphHeight+15;
			//normalize force for drawing
			min = Double.MAX_VALUE-1;
			max = -min;
			for(Boat1D bn : b){
				max = Math.max(max, force[bn.index]);
				min = Math.min(min, force[bn.index]);
			}
			if(min>=max){
				min = 0.0;
			}
			grg.setColor(Color.BLACK);
			grg.drawRect(0, posy, graphWidth-1, graphHeight);
			grg.setColor(Color.RED);
			int fpix = 0;
			for(Boat1D bn: b){
				fpix = posy +graphHeight - (int)(graphHeight*((force[bn.index]-min)/(max-min)));
				grg.fillOval((int)(bn.x/R * graphWidth), fpix-2, 5, 5);
			}
			if(min < 0 && max > 0){
				int zeroline = posy +graphHeight - (int)(graphHeight*((0-min)/(max-min)));
				grg.drawLine(0, zeroline, graphWidth, zeroline);
			}
			grg.drawString("net force on boat(g cm/s/s)", 50, posy - 2);
			grg.drawString("max = "+String.format("%.2f",max), 400, posy - 2);
			grg.drawString("min = "+String.format("%.2f",min), 600, posy - 2);
			
			//draw flow
			posy += graphHeight+15;
			//normalize flow for drawing
			min = Double.MAX_VALUE-1;
			max = -min;
			double totalfl = 0.0;
			for(Boat1D bn : b){
				max = Math.max(max, flow[bn.index]);
				min = Math.min(min, flow[bn.index]);
				totalfl += 60*bn.v;
			}
			if(min>=max){
				min = 0.0;
			}
			//max = 550; //just for convenience
			//min = 200;
			
			totalfl = totalfl/R;
			grg.setColor(Color.BLACK);
			grg.drawRect(0, posy, graphWidth-1, graphHeight);
			grg.setColor(Color.RED);
			int flpix = 0;
			for(Boat1D bn: b){
				flpix = posy +graphHeight - (int)(graphHeight*((flow[bn.index]-min)/(max-min)));
				grg.fillOval((int)(bn.x/R * graphWidth), flpix-2, 5, 5);
			}
			if(min < totalfl && max > totalfl){
				int zeroline = posy +graphHeight - (int)(graphHeight*((totalfl-min)/(max-min)));
				grg.drawLine(0, zeroline, graphWidth, zeroline);
			}
			grg.drawString("flow(boats/min) total = "+String.format("%.2f",totalfl), 50, posy - 2);
			grg.drawString("max = "+String.format("%.2f",max), 400, posy - 2);
			grg.drawString("min = "+String.format("%.2f",min), 600, posy - 2);
		}
		
		public void drawC(){
			computeData();
			double[] xval = new double[c.length];
			double[] cval = new double[c.length];
			double dx = model.getDisc()[4];
			for(int i=0; i<c.length; i++){
				xval[i] = 0.0 + i*dx;
				cval[i] = c[i];
			}
			p2d.setTitle("c(x)");
			p2d.setLabels("x", null);
			p2d.showAxes(true);
			p2d.lines(true);
			p2d.clearData();
			p2d.setXLimits(0.0, model.getDisc()[3]);
			p2d.addX(xval);
			p2d.addY(cval, p2d.DOT, Color.RED);
			p2d.paint();
		}
		
		public void drawLabelSpace(){
			if(N<1){
				return;
			}
			double[] dims = model.getParameters();
			double R = dims[0];
			double L = dims[3];
			double rc = dims[4];
			int wid = labelSpace.getWidth();
			int hig = labelSpace.getHeight();
			int segwid = (int)(wid/N);
			Graphics2D lsg = labelSpace.createGraphics();
			lsg.setColor(Color.WHITE);
			lsg.fillRect(0, 0, wid, hig);
			lsg.setColor(Color.BLACK);
			if(N < 2){
				lsg.drawString("fewer than 2 boats", 20, hig/2);
				return;
			}
			
			//draw lines for each boat
			for(int i=0; i<N; i++){
				lsg.drawLine(i*segwid, 0, i*segwid, hig);
			}
			if(velocity == null){
				return;
			}
			//draw the spacing in grayscale
			double totalspace = 0.0;
			//double maxspace = R-N*L-N*rc;
			double maxspace = 0.0;
			int ypos;
			int lineskip = velocity.length/hig;
			float rgb = 0.0f;
			for(int i=0; i<velocity.length-lineskip; i+=lineskip){
				for(int j=0; j<N; j++){
					maxspace = Math.max(maxspace, allspacing[i][j]);
				}
			}
			for(int i=0; i<velocity.length-lineskip; i+=lineskip){
				ypos = hig - (hig*i)/velocity.length;
				for(int j=0; j<allspacing[i].length; j++){
					totalspace += allspacing[i][j];
					rgb = (float)(Math.min(Math.max(allspacing[i][j]/maxspace,0),1.0));
					lsg.setColor(new Color(rgb,rgb,rgb));
					lsg.drawLine(j*segwid+1,ypos,(j+1)*segwid-1,ypos);
				}
				if(totalspace < 0.1){
					break;
				}
				totalspace = 0.0;
			}
		}
	}

	//////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	/**
     *  This is a panel for displaying and modifying parameters
     */
	class PPanel extends JPanel {
		// there are lots of labels and fields
        JLabel NLabel, RLabel, mLabel, visLabel, LLabel, rcLabel, DLabel, vDLabel, kLabel, vkLabel, aLabel, vaLabel, bLabel, csatLabel, dgLabel;
        public JTextField NField, RField, mField, visField, LField, rcField, DField, vDField, kField, vkField, aField, vaField, bField, csatField, dgField;

		public PPanel() {
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
            rcLabel = new JLabel("radius of camphor supply(cm) = ");
            rcLabel.setHorizontalAlignment(JLabel.TRAILING);
            DLabel = new JLabel("surf. diffusion(cm*cm/s) = ");
            DLabel.setHorizontalAlignment(JLabel.TRAILING);
			vDLabel = new JLabel("bulk diffusion(cm*cm/s) = ");
            vDLabel.setHorizontalAlignment(JLabel.TRAILING);
            kLabel = new JLabel("sublimation(1/s) = ");
            kLabel.setHorizontalAlignment(JLabel.TRAILING);
			vkLabel = new JLabel("surf. to vol.(1/s) = ");
            vkLabel.setHorizontalAlignment(JLabel.TRAILING);
            aLabel = new JLabel("surf. source = ");
            aLabel.setHorizontalAlignment(JLabel.TRAILING);
			vaLabel = new JLabel("bulk source = ");
            vaLabel.setHorizontalAlignment(JLabel.TRAILING);
            bLabel = new JLabel("beta = ");
            bLabel.setHorizontalAlignment(JLabel.TRAILING);
			csatLabel = new JLabel("c saturation = ");
            csatLabel.setHorizontalAlignment(JLabel.TRAILING);
			dgLabel = new JLabel("max dgamma = ");
            dgLabel.setHorizontalAlignment(JLabel.TRAILING);

			// insert the initial parameters
			double[] params =  model.getParameters();
			NField = new JTextField(String.valueOf(N), 6);
            NField.setEditable(true);
            RField = new JTextField(String.valueOf(params[0]), 6);
            RField.setEditable(true);
            mField = new JTextField(String.valueOf(params[1]), 6);
            mField.setEditable(true);
            visField = new JTextField(String.valueOf(params[2]), 6);
            visField.setEditable(true);
            LField = new JTextField(String.valueOf(params[3]), 6);
            LField.setEditable(true);
            rcField = new JTextField(String.valueOf(params[4]), 6);
            rcField.setEditable(true);
            DField = new JTextField(String.valueOf(params[5]), 6);
            DField.setEditable(true);
			vDField = new JTextField(String.valueOf(params[6]), 6);
            vDField.setEditable(true);
            kField = new JTextField(String.valueOf(params[7]), 6);
            kField.setEditable(true);
			vkField = new JTextField(String.valueOf(params[8]), 6);
            vkField.setEditable(true);
            aField = new JTextField(String.valueOf(params[9]), 6);
            aField.setEditable(true);
			vaField = new JTextField(String.valueOf(params[10]), 6);
            vaField.setEditable(true);
            bField = new JTextField(String.valueOf(params[11]), 6);
            bField.setEditable(true);
			csatField = new JTextField(String.valueOf(params[12]), 6);
            csatField.setEditable(true);
			dgField = new JTextField(String.valueOf(params[13]), 6);
            dgField.setEditable(true);

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
            this.add(rcLabel);
            this.add(rcField);
            this.add(DLabel);
            this.add(DField);
			this.add(vDLabel);
            this.add(vDField);
            this.add(kLabel);
            this.add(kField);
			this.add(vkLabel);
            this.add(vkField);
            this.add(aLabel);
            this.add(aField);
			this.add(vaLabel);
            this.add(vaField);
            this.add(bLabel);
            this.add(bField);
			this.add(csatLabel);
            this.add(csatField);
			this.add(dgLabel);
            this.add(dgField);
		}

		public void newParameters(){
			try{
				double[] params = new double[14];
				params[0] = Double.valueOf(RField.getText());
				params[1] = Double.valueOf(mField.getText());
				params[2] = Double.valueOf(visField.getText());
				params[3] = Double.valueOf(LField.getText());
				params[4] = Double.valueOf(rcField.getText());
				params[5] = Double.valueOf(DField.getText());
				params[6] = Double.valueOf(vDField.getText());
				params[7] = Double.valueOf(kField.getText());
				params[8] = Double.valueOf(vkField.getText());
				params[9] = Double.valueOf(aField.getText());
				params[10] = Double.valueOf(vaField.getText());
				params[11] = Double.valueOf(bField.getText());
				params[12] = Double.valueOf(csatField.getText());
				params[13] = Double.valueOf(dgField.getText());
				
				model.setParameters(params);
				if(model.getN() != (int)Integer.valueOf(NField.getText())){
					N = (int)Integer.valueOf(NField.getText());
					model.setN(N);
				}
			} catch(Exception expt){
				JOptionPane.showMessageDialog(this, "Couldn't read parameters. Try again.");
			}
		}
	}

}