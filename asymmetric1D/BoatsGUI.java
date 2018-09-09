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
class BoatsGUI implements ActionListener{
    // These are the parts of the GUI
	BoatsModel model;
	int N;
	ArrayList<Boat> b;
	double[] c;
	Vector<Integer> jammedBoat;
	double phasev, jamt1, jamt2, dt1, dt2;
	double[] density, dT, dc, force, flow;

	BoatsWorker worker;

	JFrame frame;
	BPanel boatsPanel;
	PPanel paramPanel;
	JPanel cntrlPanel, rightPanel, bottomPanel;
	BufferedImage spaceTime, experiment, graphs, recordst;
	Plot2D p2d;
	int imageType;
	String[] imageNames;
	int simWidth, simHeight;
	double guiTimeStart, guiTimeEnd, modelTime;

	JMenuBar menuBar;
	JMenu fileMenu, modelMenu, windowMenu, simMenu, addOnMenu;
    JMenuItem infoItem, fRateItem, fastItem, clusterItem, sim1Item;
	JMenuItem sim2Item, sim3Item, clusterTestItem, clusterTest2Item, clusterTest3Item, clusterTest4Item;
	ButtonGroup windowGroup;
	JRadioButtonMenuItem bigWindowItem, smallWindowItem, recordImageItem;

    JButton runButton, timeButton, adButton, rmButton, homoButton, randButton;
	JButton paramButton, resetButton, stopButton, displayButton, backgButton, writeButton;

	JTextField backgField;

	boolean showc;
	long compytime;
	int frameRate;
	boolean computeFast, recordNext;

	File sim1File, sim2File, imageFile;
	FileWriter writer;

	/**
	 * The constructor initializes all the model parts and prepares
	 * everything to run.
	 */
	public BoatsGUI(){
		frameRate = 30;
		showc = false;
		computeFast = false;
		recordNext = false;
		model = new BoatsModel();
		N = 0;
		b = model.getBoats();
		c = model.getc();
		jammedBoat = new Vector<Integer>();
		phasev = 0.0;
		dt1 = 0.0;
		jamt1 = 0.0;

		worker = new BoatsWorker();

		simWidth = 800;
		simHeight = 600;
		imageType = 0;
		imageNames = new String[]{"space time", "experiment", "graphs", "something"};
		
		spaceTime = new BufferedImage(simWidth, simHeight, BufferedImage.TYPE_INT_ARGB);
		spaceTime.createGraphics().setColor(Color.WHITE);
		spaceTime.createGraphics().fillRect(0, 0, simWidth, simHeight);
		
		experiment = new BufferedImage(simWidth, simHeight, BufferedImage.TYPE_INT_ARGB);
		experiment.createGraphics().setColor(Color.WHITE);
		experiment.createGraphics().fillRect(0, 0, simWidth, simHeight);
		
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

		sim1Item = new JMenuItem("vary N with fixed parameters");
		sim1Item.addActionListener(this);
		sim2Item = new JMenuItem("vary a parameter for fixed N");
		sim2Item.addActionListener(this);
		sim3Item = new JMenuItem("fix N vary param and save images");
		sim3Item.addActionListener(this);
		clusterTestItem = new JMenuItem("compute b now");
		clusterTestItem.addActionListener(this);
		clusterTest2Item = new JMenuItem("b vs parameter");
		clusterTest2Item.addActionListener(this);
		clusterTest3Item = new JMenuItem("b vs. par vs. par");
		clusterTest3Item.addActionListener(this);
		clusterTest4Item = new JMenuItem("vary N (hysteresis)");
		clusterTest4Item.addActionListener(this);

		fRateItem = new JMenuItem("frame rate");
		fRateItem.addActionListener(this);
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

		clusterItem = new JMenuItem("clusterizer");
		clusterItem.addActionListener(this);

		//put the menus together
		modelMenu.add(infoItem);
		windowMenu.add(fRateItem);
		windowMenu.add(fastItem);
		windowMenu.addSeparator();
		windowMenu.add(recordImageItem);
		windowMenu.addSeparator();
		windowMenu.add(bigWindowItem);
		windowMenu.add(smallWindowItem);
		simMenu.add(sim1Item);
		simMenu.add(sim2Item);
		simMenu.add(sim3Item);
		simMenu.add(clusterTestItem);
		simMenu.add(clusterTest2Item);
		simMenu.add(clusterTest3Item);
		simMenu.add(clusterTest4Item);
		addOnMenu.add(clusterItem);

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
			frame.repaint();
		} else if (e.getSource() == backgButton){
			model.changeBackground(Double.valueOf(backgField.getText()));
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
					writeImage(recordst, "images/N"+String.valueOf(b.size())+"D"+String.valueOf((int)pars[6])+"K"+String.valueOf((int)pars[7])+"T"+String.valueOf((int)model.getTime())+".png");
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
			computeFast = true;
			runSim1();
			computeFast = false;
		} else if (e.getSource() == sim2Item){
			computeFast = true;
			runSim2();
			computeFast = false;
		} else if (e.getSource() == sim3Item){
			runSim3();
		} else if (e.getSource() == clusterTestItem){
			System.out.println(String.valueOf(getBeta(50, false)));
		} else if (e.getSource() == clusterTest2Item){
			testForClusters();
		} else if (e.getSource() == clusterTest3Item){
			testForClustersLikeMad();
		} else if (e.getSource() == clusterTest4Item){
			testForClustersOneRun();
		} else if (e.getSource() == clusterItem){
			BoatsClusterizer clusterizer = new BoatsClusterizer(model.getParameters());
		}
	}

	/**
	 * This simulation uses the current parameters and varies N.
	 * It runs the simulation for 15 seconds for each N.
	 * It calculates: flow, ave. c, delta c, ave. v, delta v, number of jams
	 * And writes them to a file "data/sim1" in columns
	 * Also writes a shell file to run gnuplot
	 */
	public void runSim1(){
		//we need several numbers
		double[] flow, cave, dc, vave, dv, jnum;
		double R = model.getDisc()[3];
		int nx = (int)model.getDisc()[5];
		int maxN = Integer.valueOf((String)JOptionPane.showInputDialog(frame, "Input maximum number of boats"));
		double maxT = Double.valueOf((String)JOptionPane.showInputDialog(frame, "Input time interval"));

		flow = new double[maxN];
		cave = new double[maxN];
		dc = new double[maxN];
		vave = new double[maxN];
		dv = new double[maxN];
		jnum = new double[maxN];

		// set T=15 sec.
		model.setInterval(maxT);
		guiTimeStart = model.getTime();
		guiTimeEnd = guiTimeStart+maxT;

		//loop over each N
		for(int i=0; i<maxN; i++){
			reset();

			N = i+1;
			model.setN(N);
			paramPanel.NField.setText(String.valueOf(N));
			model.randomize();

			worker = new BoatsWorker();
			worker.execute();
			while(!worker.isDone()){
				//boatsPanel.paintImmediately(0,0,boatsPanel.getWidth(),boatsPanel.getHeight());
				try{
					synchronized(this){
						wait(100);
					}
				}catch(InterruptedException ignore){}
			}
			System.out.println("finished N = "+String.valueOf(N));
			boatsPanel.paintImmediately(0,0,boatsPanel.getWidth(),boatsPanel.getHeight());

			// compute flow, vave and dv
			flow[i] = 0.0;
			vave[i] = 0.0;
			dv[i] = 0.0;
			double maxv = 0.0;
			double minv = 10.0;
			for (int j=0; j<N; j++){
				vave[i] += b.get(j).v;
				maxv = Math.max(maxv, b.get(j).v);
				minv = Math.min(minv, b.get(j).v);
			}
			vave[i] = vave[i]/N;
			flow[i] = vave[i]*N/R*60;
			dv[i] = maxv-minv;

			// compute cave and dc
			cave[i] = 0.0;
			dc[i] = 0.0;
			double maxc = 0.0;
			double minc = 100.0;
			for (int j=0; j<nx; j++){
				cave[i] += c[j];
				maxc = Math.max(maxc, c[j]);
				minc = Math.min(minc, c[j]);
			}
			cave[i] = cave[i]/nx;
			dc[i] = maxc-minc;

			//compute jnum (number of jams perhaps)
			jnum[i] = 0;
			if(dv[i] > vave[i]*0.3){
				boolean up = (b.get(0).v > vave[i]+dv[i]/4);
				boolean change = !up;
				for (int j=0; j<N; j++){
					if(b.get(j).v > vave[i]+dv[i]/4){
						up = true;
						change = false;
					}
					if(b.get(j).v < vave[i]-dv[i]/4){
						up = false;
					}
					if(up == change){
						change = true;
						jnum[i]++;
					}
				}
				
				jnum[i] = jammedBoat.size();
			}
		}

		//now write it to a file in columns
		sim1File = new File("data");
		if(!sim1File.exists()){
			sim1File.mkdir();
		}
		//sim1File = new File("data\\sim1");	//for windows
		sim1File = new File("data/sim1");	//for everything else
		try{
			writer = new FileWriter(sim1File);
			for(int i=0; i<maxN; i++){
				writer.write(String.valueOf(i+1)+" "+String.valueOf(flow[i])+" "+String.valueOf(vave[i])+" "+String.valueOf(dv[i])+" "+String.valueOf(cave[i])+" "+String.valueOf(dc[i])+" "+String.valueOf(jnum[i])+"\n");
			}
			JOptionPane.showMessageDialog(frame, "successfully recorded data");
			writer.close();
		}catch(Exception e){
			System.err.println("Error: " + e.getMessage());
			JOptionPane.showMessageDialog(frame, "file writing error");
		}
	}
	
	
	/**
	 * This keeps N constant and varies some parameter
	 */
	public void runSim2(){
		if(JOptionPane.showConfirmDialog(frame, "First set the parameters to their lowest values\nAre they set?", "choose one", JOptionPane.YES_NO_OPTION) == JOptionPane.NO_OPTION){
			return;
		}
		Object[] choices = {"K", "alpha", "beta"};
		Object choice = JOptionPane.showInputDialog(frame, "Choose the parameter to vary:", "setup simulation", JOptionPane.PLAIN_MESSAGE, null, choices, "ham");
		if(choice == null){
			return;
		}
		double maxp = Double.valueOf((String)JOptionPane.showInputDialog(frame, "Input maximum value of the parameter"));
		double maxT = Double.valueOf((String)JOptionPane.showInputDialog(frame, "Input time interval to run each simulation(seconds)"));

		//we need several numbers
		double[] flow, cave, dc, vave, dv, jnum;
		double R = model.getDisc()[3];
		int nx = (int)model.getDisc()[5];
		int points = 20;
		double[] para = new double[points];
		if(choice == choices[0]){
			para[0] = model.getParameters()[7];
		}else if(choice == choices[1]){
			para[0] = model.getParameters()[8];
		}else if(choice == choices[2]){
			para[0] = model.getParameters()[9];
		}
		double dp = (maxp - para[0])/(points-1);
		for(int i=1; i<points; i++){
			para[i] = para[i-1]+dp;
		}

		flow = new double[points];
		cave = new double[points];
		dc = new double[points];
		vave = new double[points];
		dv = new double[points];
		jnum = new double[points];

		// set T
		model.setInterval(maxT);
		guiTimeStart = model.getTime();
		guiTimeEnd = guiTimeStart+maxT;

		//loop over each data point
		for(int i=0; i<points; i++){
			reset();
			if(choice == choices[0]){
				paramPanel.kField.setText(String.valueOf(para[i]));
				paramPanel.newParameters();
			}else if(choice == choices[1]){
				paramPanel.aField.setText(String.valueOf(para[i]));
				paramPanel.newParameters();
			}else if(choice == choices[2]){
				paramPanel.bField.setText(String.valueOf(para[i]));
				paramPanel.newParameters();
			}
			model.randomize();

			worker = new BoatsWorker();
			worker.execute();
			while(!worker.isDone()){
				//boatsPanel.paintImmediately(0,0,boatsPanel.getWidth(),boatsPanel.getHeight());
				try{
					synchronized(this){
						wait(100);
					}
				}catch(InterruptedException ignore){}
			}
			System.out.println("finished data point "+String.valueOf(i));

			// compute flow, vave and dv
			flow[i] = 0.0;
			vave[i] = 0.0;
			dv[i] = 0.0;
			double maxv = 0.0;
			double minv = 10.0;
			for (int j=0; j<N; j++){
				vave[i] += b.get(j).v;
				maxv = Math.max(maxv, b.get(j).v);
				minv = Math.min(minv, b.get(j).v);
			}
			vave[i] = vave[i]/N;
			flow[i] = vave[i]*N/R*60;
			dv[i] = maxv-minv;

			// compute cave and dc
			cave[i] = 0.0;
			dc[i] = 0.0;
			double maxc = 0.0;
			double minc = 100.0;
			for (int j=0; j<nx; j++){
				cave[i] += c[j];
				maxc = Math.max(maxc, c[j]);
				minc = Math.min(minc, c[j]);
			}
			cave[i] = cave[i]/nx;
			dc[i] = maxc-minc;

			//compute jnum (number of jams perhaps)
			jnum[i] = 0;
			if(dv[i] > vave[i]*0.3){
				boolean up = (b.get(0).v > vave[i]+dv[i]/4);
				boolean change = !up;
				for (int j=0; j<N; j++){
					if(b.get(j).v > vave[i]+dv[i]/4){
						up = true;
						change = false;
					}
					if(b.get(j).v < vave[i]-dv[i]/4){
						up = false;
					}
					if(up == change){
						change = true;
						jnum[i]++;
					}
				}
				
				jnum[i] = jammedBoat.size();
			}
		}

		//now write it to a file in columns
		sim2File = new File("data");
		if(!sim2File.exists()){
			sim2File.mkdir();
		}
		//sim1File = new File("data\\sim2"+(String)choice);	//for windows
		sim2File = new File("data/sim2"+(String)choice);	//for everything else
		try{
			writer = new FileWriter(sim2File);
			for(int i=0; i<points; i++){
				writer.write(String.valueOf(para[i])+" "+String.valueOf(flow[i])+" "+String.valueOf(vave[i])+" "+String.valueOf(dv[i])+" "+String.valueOf(cave[i])+" "+String.valueOf(dc[i])+" "+String.valueOf(jnum[i])+"\n");
			}
			JOptionPane.showMessageDialog(frame, "successfully recorded data");
			writer.close();
		}catch(Exception e){
			System.err.println("Error: " + e.getMessage());
			JOptionPane.showMessageDialog(frame, "file writing error");
		}
	}
	
	/**
	 * This keeps N constant and varies some parameter and saves images
	 */
	public void runSim3(){
		if(JOptionPane.showConfirmDialog(frame, "First set the parameters to their lowest values\nAre they set?", "choose one", JOptionPane.YES_NO_OPTION) == JOptionPane.NO_OPTION){
			return;
		}
		Object[] choices = {"D", "K", "beta"};
		Object choice = JOptionPane.showInputDialog(frame, "Choose the parameter to vary:", "setup simulation", JOptionPane.PLAIN_MESSAGE, null, choices, "ham");
		if(choice == null){
			return;
		}
		double maxp = Double.valueOf((String)JOptionPane.showInputDialog(frame, "Input maximum value of the parameter"));
		int imageCount = Integer.valueOf((String)JOptionPane.showInputDialog(frame, "Input number of images to save"));
		double maxT = Double.valueOf((String)JOptionPane.showInputDialog(frame, "Input time interval to run each simulation(seconds)"));
		
		//we need several numbers
		int number = model.getN();
		double[] flow, cave, dc, vave, dv, jnum;
		int points = imageCount;
		double[] para = new double[points];
		double R = model.getDisc()[3];
		int nx = (int)model.getDisc()[5];
		if(choice == choices[0]){
			para[0] = model.getParameters()[6];
		}else if(choice == choices[1]){
			para[0] = model.getParameters()[7];
		}else if(choice == choices[2]){
			para[0] = model.getParameters()[9];
		}
		double dp = (maxp - para[0])/(points-1);
		for(int i=1; i<points; i++){
			para[i] = para[i-1]+dp;
		}
		
		flow = new double[points];
		cave = new double[points];
		dc = new double[points];
		vave = new double[points];
		dv = new double[points];
		jnum = new double[points];
		
		// set T
		model.setInterval(maxT);
		
		recordNext = true;
		frameRate = 1;
		
		//loop over each data point
		for(int i=0; i<points; i++){
			reset();
			
			if(choice == choices[0]){
				paramPanel.DField.setText(String.valueOf(para[i]));
				paramPanel.NField.setText(String.valueOf(number));
				paramPanel.newParameters();
			}else if(choice == choices[1]){
				paramPanel.kField.setText(String.valueOf(para[i]));
				paramPanel.NField.setText(String.valueOf(number));
				paramPanel.newParameters();
			}else if(choice == choices[2]){
				paramPanel.bField.setText(String.valueOf(para[i]));
				paramPanel.NField.setText(String.valueOf(number));
				paramPanel.newParameters();
			}
			
			model.randomize();
			if(worker.isDone()){
				worker = new BoatsWorker();
			}
			worker.execute();
			while(!worker.isDone()){
				//boatsPanel.paintImmediately(0,0,boatsPanel.getWidth(),boatsPanel.getHeight());
				try{
					Thread.currentThread().sleep(100);
				}catch(InterruptedException ignore){}
			}
			System.out.println("finished data point "+String.valueOf(i));
			worker.done();
			
			b = model.getBoats();
			c = model.getc();
			
			// compute flow, vave and dv
			flow[i] = 0.0;
			vave[i] = 0.0;
			dv[i] = 0.0;
			double maxv = 0.0;
			double minv = 10.0;
			for (int j=0; j<N; j++){
				vave[i] += b.get(j).v;
				maxv = Math.max(maxv, b.get(j).v);
				minv = Math.min(minv, b.get(j).v);
			}
			vave[i] = vave[i]/N;
			flow[i] = vave[i]*N/R*60;
			dv[i] = maxv-minv;
			
			// compute cave and dc
			cave[i] = 0.0;
			dc[i] = 0.0;
			double maxc = 0.0;
			double minc = 100.0;
			for (int j=0; j<nx; j++){
				cave[i] += c[j];
				maxc = Math.max(maxc, c[j]);
				minc = Math.min(minc, c[j]);
			}
			cave[i] = cave[i]/nx;
			dc[i] = maxc-minc;
			
			//compute jnum (number of jams perhaps)
			jnum[i] = 0;
			if(dv[i] > vave[i]*0.3){
				boolean up = (b.get(0).v > vave[i]+dv[i]/4);
				boolean change = !up;
				for (int j=0; j<N; j++){
					if(b.get(j).v > vave[i]+dv[i]/4){
						up = true;
						change = false;
					}
					if(b.get(j).v < vave[i]-dv[i]/4){
						up = false;
					}
					if(up == change){
						change = true;
						jnum[i]++;
					}
				}
				
				jnum[i] = jammedBoat.size();
			}
			
			//take a breath and let the images finish drawing
//			try{
//				Thread.currentThread().sleep(1500);
//			}catch(InterruptedException ignore){}
			
			double[] pars = model. getParameters();
			writeImage(recordst, "images/N"+String.valueOf(b.size())+"D"+String.valueOf((int)(10*pars[6]))+"K"+String.valueOf((int)pars[7])+"T"+String.valueOf((int)model.getTime())+String.valueOf(System.currentTimeMillis())+".png");
			recordst.flush();
			recordst = null;
		}
		
		recordNext = false;
		frameRate = 30;
		
		//now write it to a file in columns
		double[] pars = model. getParameters();
		sim2File = new File("data");
		if(!sim2File.exists()){
			sim2File.mkdir();
		}
		sim2File = new File("data/sim2"+"N"+String.valueOf(b.size())+"D"+String.valueOf((int)pars[6])+"K"+String.valueOf((int)pars[7])+"T"+String.valueOf((int)model.getTime()));	//for everything else
		try{
			writer = new FileWriter(sim2File);
			for(int i=0; i<points; i++){
				writer.write(String.valueOf(para[i])+" "+String.valueOf(flow[i])+" "+String.valueOf(vave[i])+" "+String.valueOf(dv[i])+" "+String.valueOf(cave[i])+" "+String.valueOf(dc[i])+" "+String.valueOf(jnum[i])+"\n");
			}
			JOptionPane.showMessageDialog(frame, "successfully recorded data");
			writer.close();
		}catch(Exception e){
			System.err.println("Error: " + e.getMessage());
			JOptionPane.showMessageDialog(frame, "file writing error");
		}
	}
	
	public void testForClusters(){
		// run several tests and display the data
		int points = 15;
		double[] K = new double[points];
		for(int i=0; i<points; i++){
			K[i] = i*0.05 + model.getParameters()[7];
		}
		double[] beta = new double[points];
		for(int i=0; i<points; i++){
			reset();
			paramPanel.NField.setText(String.valueOf(10));
			paramPanel.kField.setText(String.valueOf(K[i]));
			paramButton.doClick();
			model.randomize();
			model.setInterval(60);
			guiTimeStart = model.getTime();
			guiTimeEnd = guiTimeStart+60;
			computeFast = true;
			
			if(worker.isDone()){
				worker = new BoatsWorker();
			}
			worker.execute();
			while(!worker.isDone()){
				//boatsPanel.paintImmediately(0,0,boatsPanel.getWidth(),boatsPanel.getHeight());
				try{
					Thread.currentThread().sleep(100);
				}catch(InterruptedException ignore){}
			}
			
			beta[i] = getBeta(20, false);
			
		}
		
		Plot2D betaDplot = new Plot2D(700, 500);
		betaDplot.showAxes(true);
		betaDplot.setLabels("K", "<beta>D");
		betaDplot.boxed(true);
		betaDplot.lines(true);
		betaDplot.addX(K);
		betaDplot.addY(beta);
		betaDplot.setXLimits(K[0], K[K.length-1]);
		
		betaDplot.paint();
		
		PlotWindow pwind = new PlotWindow(betaDplot);
	}
	
	/*
	 / compute beta for varying N (or one parameter)
	 / Test for hysteresis. write data to file and draw plot
	 */
	public void testForClustersOneRun(){
		// run one long simulation with changing stuffs and record at each change
		int points = 20;
		
		double relaxT = 100;
		double recordT = 100;
		boolean goup = JOptionPane.showConfirmDialog(frame, "choose yes for up, no for down", "choose one", JOptionPane.YES_NO_OPTION) == JOptionPane.YES_OPTION;
		try{
			points = Integer.valueOf((String)JOptionPane.showInputDialog(frame, "How many data points?"));
			relaxT = Double.valueOf((String)JOptionPane.showInputDialog(frame, "Input time interval to run each simulation between recording"));
			recordT = Double.valueOf((String)JOptionPane.showInputDialog(frame, "Input time interval to record"));
		}catch(Exception e){
			return;
		}
		double[] beta = new double[points];
		double[] par = new double[points];
		int tempN = N;
		double tempk = model.getParameters()[7];
		double tempa = model.getParameters()[8];
		double dk = 0.025;
		double da = 0.2;
		
		int variable = 2; //0=N, 1=k, 2=alpha
		
		computeFast = true;
		// set up an initial state
		model.randomize();
		// do all the time steps the hard way
		int stepCount = (int)model.getDisc()[2]; 
		for(int j=0; j<stepCount; j++){
			model.timeStep();
		}
		
		if(variable == 0){
			// run for each value of N
			if(goup){ rmButton.doClick(); }
			if(!goup){ adButton.doClick(); }
		}else if(variable == 1){
			for(int i=0; i<points; i++){
				if(goup){
					par[i] = tempk + i*dk;
				}else{
					par[i] = tempk - i*dk;
				}
			}
		}else if(variable == 2){
			for(int i=0; i<points; i++){
				if(goup){
					par[i] = tempa + i*da;
				}else{
					par[i] = tempa - i*da;
				}
			}
		}
		
		for(int i=0; i<points; i++){
			if(variable == 0){
				if(goup){
					adButton.doClick();
				}else{
					rmButton.doClick();
					homoButton.doClick();
				}
				par[i] = N;
			}else if(variable == 1){
				paramPanel.kField.setText(String.valueOf(par[i]));
				paramPanel.newParameters();
			}else if(variable == 2){
				paramPanel.aField.setText(String.valueOf(par[i]));
				paramPanel.newParameters();
			}
			
			
			// first relax the system
			model.setInterval(relaxT);
			guiTimeStart = model.getTime();
			guiTimeEnd = guiTimeStart+relaxT;
			
			// do all the time steps the hard way
			stepCount = (int)model.getDisc()[2]; 
			for(int j=0; j<stepCount; j++){
				model.timeStep();
			}
			/*
			temptime = model.getTime();
			while(model.getTime() < temptime+relaxT-10){
				model.setInterval(relaxT);
				guiTimeStart = model.getTime();
				guiTimeEnd = guiTimeStart+relaxT;
				
				worker = new BoatsWorker();
				worker.execute();
				while(!worker.isDone()){
					try{
						synchronized(this){
							wait(100);
						}
					}catch(InterruptedException ignore){
						System.out.println("interrupted");
					}
				}
			}
			 */
			
			// then record data
			beta[i] = getBeta(recordT, false);
			
			System.out.println("Finished "+String.valueOf(i+1)+" of "+String.valueOf(points)+" b="+String.valueOf(beta[i]));
		}
		
		// plot it
		Plot2D betaDplot = new Plot2D(700, 500);
		betaDplot.showAxes(true);
		betaDplot.setLabels("par", "b");
		betaDplot.boxed(true);
		betaDplot.lines(true);
		betaDplot.addX(par);
		betaDplot.addY(beta);
		betaDplot.setXLimits(par[0], par[par.length-1]);
		betaDplot.paint();
		PlotWindow pwind = new PlotWindow(betaDplot);
		
		// record it
		sim2File = new File("hysteresis");
		if(!sim2File.exists()){
			sim2File.mkdir();
		}
		sim2File = new File("hysteresis/hysN"+String.valueOf(tempN)+"K"+String.valueOf((int)tempk)+"t"+String.valueOf(System.currentTimeMillis()));
		try{
			writer = new FileWriter(sim2File);
			writer.write("# par, b\n");
			for(int i=0; i<points; i++){
				writer.write(String.valueOf(par[i])+" "+String.valueOf(beta[i])+"\n");
			}
			writer.close();
		}catch(Exception e){
			System.err.println("Error: " + e.getMessage());
			JOptionPane.showMessageDialog(frame, "file writing error");
		}
		
		//computeFast = false;
	}
	
	/* 
	 / compute b for many things
	 / redord data
	 / and stuff
	 */
	public void testForClustersLikeMad(){
		// run several tests and display the data
		long startingTime = System.currentTimeMillis();
		int points = 20;
		double[] K = new double[points];
		double[] M = new double[points];
		double[] G = new double[points];
		double[] alph = new double[points];
		double[] D = new double[points];
		int whichParam = 4;	//0=vis, 1=gamma, 2=alpha, 3=D, 4 = N
		for(int i=0; i<points; i++){
			K[i] = i*0.1 + model.getParameters()[7];
			M[i] = i*0.003 + model.getParameters()[2];
			G[i] = i*0.05 - 0.3 + model.getParameters()[3];
			alph[i] = i*5 + model.getParameters()[8];
			D[i] = i*.5 + model.getParameters()[6];
		}
		double[] beta = new double[points];
		
		//write it to a file in matrix form
		sim2File = new File("clusterdata");
		if(!sim2File.exists()){
			sim2File.mkdir();
		}
		if(whichParam == 0){
			sim2File = new File("clusterdata/varyKandM");
			try{
				writer = new FileWriter(sim2File);
				writer.write("# K="+String.valueOf(K[0])+" to "+String.valueOf(K[K.length-1])+" , M="+String.valueOf(M[0])+" to "+String.valueOf(M[M.length-1])+"\n");
				writer.write("# K, M, b\n\n");
			}catch(Exception e){
				System.err.println("Error: " + e.getMessage());
				JOptionPane.showMessageDialog(frame, "file writing error");
			}
		}else if(whichParam == 1){
			sim2File = new File("clusterdata/varyKandG");
			try{
				writer = new FileWriter(sim2File);
				writer.write("# K="+String.valueOf(K[0])+" to "+String.valueOf(K[K.length-1])+" , G="+String.valueOf(G[0])+" to "+String.valueOf(G[G.length-1])+"\n");
				writer.write("# K, G, b\n\n");
			}catch(Exception e){
				System.err.println("Error: " + e.getMessage());
				JOptionPane.showMessageDialog(frame, "file writing error");
			}
		}else if(whichParam == 2){
			sim2File = new File("clusterdata/varyKandalpha");
			try{
				writer = new FileWriter(sim2File);
				writer.write("# K="+String.valueOf(K[0])+" to "+String.valueOf(K[K.length-1])+" , alpha="+String.valueOf(alph[0])+" to "+String.valueOf(alph[alph.length-1])+"\n");
				writer.write("# K, alpha, b\n\n");
			}catch(Exception e){
				System.err.println("Error: " + e.getMessage());
				JOptionPane.showMessageDialog(frame, "file writing error");
			}
		}else if(whichParam == 3){
			sim2File = new File("clusterdata/varyKandD");
			try{
				writer = new FileWriter(sim2File);
				writer.write("# K="+String.valueOf(K[0])+" to "+String.valueOf(K[K.length-1])+" , D="+String.valueOf(D[0])+" to "+String.valueOf(D[D.length-1])+"\n");
				writer.write("# K, D, b\n\n");
			}catch(Exception e){
				System.err.println("Error: " + e.getMessage());
				JOptionPane.showMessageDialog(frame, "file writing error");
			}
		}else if(whichParam == 4){
			sim2File = new File("clusterdata/varyKandN");
			try{
				writer = new FileWriter(sim2File);
				writer.write("# K="+String.valueOf(K[0])+" to "+String.valueOf(K[K.length-1])+" , N="+String.valueOf(5)+" to "+String.valueOf(45)+"\n");
				writer.write("# K, N, b\n\n");
			}catch(Exception e){
				System.err.println("Error: " + e.getMessage());
				JOptionPane.showMessageDialog(frame, "file writing error");
			}
		}
		
		for(int i=0; i<points; i++){
			for(int j=0; j<points; j++){
				reset();
				paramPanel.NField.setText(String.valueOf(10));
				paramPanel.kField.setText(String.valueOf(K[j]));
				if(whichParam == 2){
					paramPanel.aField.setText(String.valueOf(alph[i]));
				}else if(whichParam == 0){
					paramPanel.visField.setText(String.valueOf(M[i]));
				}else if(whichParam == 1){
					//model.setl(G[i]);
					System.err.println("Error: oops, remember you killed this");
				}else if(whichParam == 3){
					paramPanel.DField.setText(String.valueOf(D[i]));
				}else if(whichParam == 4){
					paramPanel.NField.setText(String.valueOf(5+i*2));
				}
				paramButton.doClick();
				model.randomize();
				model.setInterval(60);
				guiTimeStart = model.getTime();
				guiTimeEnd = guiTimeStart+60;
				computeFast = true;
				
				if(worker.isDone()){
					worker = new BoatsWorker();
				}
				worker.execute();
				while(!worker.isDone()){
					//boatsPanel.paintImmediately(0,0,boatsPanel.getWidth(),boatsPanel.getHeight());
					try{
						Thread.currentThread().sleep(100);
					}catch(InterruptedException ignore){}
				}
				
				beta[j] = getBeta(15, false);
				
				try{
					writer.write(String.valueOf(K[j])+" ");
					if(whichParam == 2){
						writer.write(String.valueOf(alph[i])+" ");
					}else if(whichParam == 0){
						writer.write(String.valueOf(M[i])+" ");
					}else if(whichParam == 1){
						writer.write(String.valueOf(G[i]*22/0.6)+" ");
					}else if(whichParam == 3){
						writer.write(String.valueOf(D[i])+" ");
					}else if(whichParam == 4){
						writer.write(String.valueOf(N)+" ");
					}
					writer.write(String.valueOf(Math.max(beta[j],0.0))+"\n");
				}catch(Exception e){
					System.err.println("Error: " + e.getMessage());
				}
			}
			try{
				writer.write("\n");
			}catch(Exception e){
				System.err.println("Error: " + e.getMessage());
			}
		}
		try{
			writer.close();
		}catch(Exception e){
			System.err.println("Error: " + e.getMessage());
		}
		long totalTime = System.currentTimeMillis() - startingTime;
		System.out.println("finished! total time = "+String.valueOf(totalTime/1000)+" seconds");
	}

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
		
		int indb, indf;
		double Tb, Tf;
		density = new double[b.size()];
		dT = new double[b.size()];
		dc = new double[b.size()];
		force = new double[b.size()];
		flow = new double[b.size()];
		
		for(Boat bn: b){
			indb = (int)Math.round(bn.x/dx);		// find the index for the back
			indf = (int)Math.round((bn.x+L)%R/dx);		// find the index for the front
			Tb = 22/((beta*c[indb])*(beta*c[indb]) + 1) + 50; //surface tension at front
			Tf = 22/((beta*c[indf])*(beta*c[indf]) + 1) + 50; //surface tension at back
			
			density[bn.index] = 2/(bn.back+bn.front+2*(L+rc));
			dT[bn.index] = Tf-Tb;
			dc[bn.index] = c[indb] - c[indf];
			force[bn.index] = (Tf-Tb)*l - vis*bn.v;
			flow[bn.index] = 60*bn.v*density[bn.index];
		}
	}
	
	public double getBeta(double interval, boolean doplot){
		//Note: this does not use the worker thread
		//the gui will be frozen like a biotch
		model.setInterval(interval);
		guiTimeStart = model.getTime();
		guiTimeEnd = guiTimeStart+interval;
		int stepCount = (int)model.getDisc()[2];
		double dt = model.getDisc()[1];
		
		BoatsClusterAnalyzer bca = new BoatsClusterAnalyzer(N, (int)(stepCount/4), 4.0*dt, 1.0/(model.getDisc()[3]-N*model.getParameters()[4]));
		double[] tempx = new double[N];
		double[] tempfront = new double[N];
		// do all the time steps and save data
		for(int i=0; i<stepCount; i++){
			model.timeStep();
			if(i%4 == 0){
				b = model.getBoats();
				for(int j=0; j<N; j++){
					tempx[j] = b.get(j).x;
					tempfront[j] = b.get(j).front;
				}
				bca.addData(tempx, tempfront);
			}
		}
		
		if(doplot){
			bca.plotBetaOfT();
		}
		
		return bca.computeMaxBeta();
	}
	
	public void writeData(){
		sim1File = new File("data");
		if(!sim1File.exists()){
			sim1File.mkdir();
		}
		//sim1File = new File("data\\captureN"+String.valueOf(boats.size()));	//for windows
		sim1File = new File("data/capture"+String.valueOf(b.size())+String.valueOf((int)phasev));	//for everything else
		try{
			writer = new FileWriter(sim1File);
			for(Boat bn: b){
				writer.write(String.valueOf(bn.x)+" "+String.valueOf(dT[bn.index])+" "+String.valueOf(dc[bn.index])+" "+String.valueOf(force[bn.index])+" "+String.valueOf(flow[bn.index])+"\n");
			}
			JOptionPane.showMessageDialog(frame, "successfully recorded data\nFile: "+sim1File.getName()+"\ncolumns: x, density, dT, dc, force, flow");
			writer.close();
		}catch(Exception e){
			System.err.println("Error: " + e.getMessage());
			JOptionPane.showMessageDialog(frame, "file writing error");
		}
	}
	
	public void writeImage(BufferedImage img, String name){
		if(img != null){
			imageFile = new File("images");
			if(!imageFile.exists()){
				imageFile.mkdir();
			}
			double[] pars = model. getParameters();
			imageFile = new File(name);	//for everything else
			try{
				ImageIO.write(img, "png", imageFile);
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
		p2d= new Plot2D(simWidth, simHeight);
		p2d.setTitle("c(x)");
		p2d.setLabels("x", null);
		p2d.setXLimits(0.0, 45.5);
		p2d.showAxes(true);
		p2d.lines(true);

		jammedBoat.clear();

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
		Thread thread;
		boolean drawn;
		BufferedImage recordtemp;

		@Override
        protected String doInBackground() {
			thread = Thread.currentThread();
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

			ArrayList<Boat> bgb = new ArrayList<Boat>(N);
			double[] bgc = new double[0];
			
			//for an extra fast computation without the fancy stuff
			if(computeFast){
				// do all the time steps and nothing else
				for(int i=0; i<disc[2]; i++){
					if(isCancelled()){
						return "cancelled";
					}
					model.timeStep();
				}
				stTime += totalFrames*disc[1];
				//bgst.createGraphics().drawString("space time diagram not updated in fast compute mode", 10, 200);
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
					
					// this colors c on the st diagram
					if(showc){
						//color the stdiagram with the c concentration
						//normalize c for drawing
						for(Double cx : bgc){
							if(cx>maxc)
								maxc = cx.floatValue();
							if(cx<minc)
								minc = cx.floatValue();
						}
						if(minc>=maxc){
							maxc = 1;
							minc = 0;
						}
						for(int k=0; k<disc[5]; k++){
							cColor[k] = Color.getHSBColor((float)(0.75-0.75*((bgc[k]-minc)/(maxc-minc))), (float)0.3, (float)1.0);
							if(bgst.getRGB(2*k, bgst.getHeight()-tpix) != Color.RED.getRGB()){
								bgst.setRGB(2*k, bgst.getHeight()-tpix, cColor[k].getRGB());
							}
							if(bgst.getRGB(2*k+1, bgst.getHeight()-tpix) != Color.RED.getRGB()){
								bgst.setRGB(2*k+1, bgst.getHeight()-tpix, cColor[k].getRGB());
							}
						}
					}
					
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
			b = (ArrayList<Boat>) billy[0];
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
			
			jammedBoat.clear();
			double maxv = -100.0;
			double minv = 100.0;
			for(Boat bn: b){
				maxv = Math.max(maxv, bn.v);
				minv = Math.min(minv, bn.v);
			}
			for(int i=1; i<b.size(); i++){
				if(b.get(i-1).v - b.get(i).v > 0.3*(maxv-minv)){
					jammedBoat.add(new Integer(i));
				}
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
					try{
						Thread.currentThread().sleep(500);
					}catch(InterruptedException ignore){}
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
			double l = dims[3];
			double L = dims[4];
			double rc = dims[5];

			int esimHeight = (int)(l/R*simWidth) + 2;
			int simposx = 20;
			int simposy = height-70;

			int picposx = simposx;
			int picposy = 20;

			//normalize c for drawing
			float minc = 100;
			float maxc = 0;
			for(Double cx : c){
				if(cx>maxc)
					maxc = cx.floatValue();
				if(cx<minc)
					minc = cx.floatValue();
			}
			if(minc>=maxc){
				minc = 0;
			}
			//normalize v for drawing
			float minv = 100;
			float maxv = 0;
			for(Boat bn : b){
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
			ig.setColor(Color.RED);
			int cpix = 0;
			for(int i=0; i<c.length; i++){
				cpix = simposy - 10 - (int)(100*((c[i]-minc)/(maxc-minc)));
				ig.drawLine((int)(simposx+i*800/c.length), cpix, (int)(simposx+i*800/c.length), cpix);
			}
			ig.drawString("camphor concentration, threshold = "+String.format("%.2f",1/dims[9]), simposx+100, simposy-115);
			ig.drawString(String.format("%.2f",maxc), simWidth+simposx+5, simposy-100);
			ig.drawString(String.format("%.2f",minc), simWidth+simposx+5, simposy-10);

			//draw the velocity
			ig.setColor(Color.BLUE);
			int vpix = 0;
			for(Boat bn: b){
				vpix = simposy - 10 - (int)(100*((bn.v-minv)/(maxv-minv)));
				ig.fillOval(simposx+(int)(bn.x/R * simWidth), vpix-2, 5, 5);
			}
			ig.drawString("velocity(cm/s)", simposx+400, simposy-115);
			ig.drawString(String.format("%.2f",maxv), simWidth+simposx+45, simposy-100);
			ig.drawString(String.format("%.2f",minv), simWidth+simposx+45, simposy-10);

			//draw the boats
			ig.setColor(Color.BLACK);
			int csize = (int)(2*rc/R*simWidth);
			int bsize = (int)(L/R*simWidth);
			for(Boat bn: b){
				int bpos = (int) (bn.x/R * simWidth);
				int cpos = (int) ((bn.x-rc)/R * simWidth);
				ig.fillOval(bpos+simposx, simposy+(int)(esimHeight/2 - bsize/2), bsize, bsize);
				ig.drawOval(cpos+simposx, simposy+(int)(esimHeight/2 - csize/2), csize, csize);
				if(bn.index == 0){
					ig.setColor(Color.RED);
					ig.fillOval(bpos+simposx, simposy+(int)(esimHeight/2 - bsize/2), bsize, bsize);
					ig.setColor(Color.BLACK);
				}
			}
			ig.setColor(Color. BLUE);
			for(Integer jn: jammedBoat){
				if(jn.intValue() < b.size()){
					int bpos = (int) (b.get(jn.intValue()).x/R * simWidth);
					ig.fillOval(bpos+simposx, simposy+(int)(esimHeight/2 - bsize/2), bsize, bsize);
				}
			}
			
			//draw the phase velocity
			dt2 = modelTime;
			if(dt2-dt1 > 0.4 && jammedBoat.size() > 0){
				double framedt = dt2-dt1;
				jamt2 = b.get(jammedBoat.get(0).intValue()).x;
				phasev = (jamt2-jamt1)/(dt2-dt1);
				jamt1 = jamt2;
				dt1 = dt2;
			}
			
			ig.setColor(Color.BLACK);
			ig.drawString("phase velocity for one jam(cm/s) = "+String.format("%.2f",phasev), simposx+400, simposy-127);

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
				computeData();
				double[] xval = new double[c.length];
				double[] cval = new double[c.length];
				double dx = model.getDisc()[4];
				int i = 0;
				for(Double cn: c){
					xval[i] = 0.0 + i*dx;
					cval[i] = cn.doubleValue();
					i++;
				}
				p2d.clearData();
				p2d.setXLimits(0.0, model.getDisc()[3]);
				p2d.addX(xval);
				p2d.addY(cval, p2d.DOT, Color.RED);
				p2d.paint();
				ig.drawImage(p2d, picposx, picposy, this);
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
			double l = dims[3];
			double L = dims[4];
			double rc = dims[5];
			Graphics2D exg = experiment.createGraphics();
			int exDia = experiment.getHeight()-40;
			double scale = (exDia)/(R/Math.PI + l);
			int boatDia = (int)(l*scale);
			exg.setColor(Color.WHITE);
			exg.fillRect(0, 0, experiment.getWidth(), experiment.getHeight());
			exg.setColor(Color.BLACK);
			exg.drawOval(20, 20, exDia, exDia);
			exg.drawOval(20+2*boatDia, 20+2*boatDia, exDia-4*boatDia, exDia-4*boatDia);
			exg.drawLine(50+exDia, 20, 50+exDia, 20+exDia);
			exg.drawLine(40+exDia, 20, 60+exDia, 20);
			exg.drawLine(40+exDia, 20+exDia, 60+exDia, 20+exDia);
			exg.drawString(String.format("%.2f", R/Math.PI)+"cm", 55+exDia, 10+exDia/2);
			for(Boat bn: b){
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
			for(Boat bn: b){
				xval[bn.index] = bn.x;
			}

			computeData();
			
			//clear the image
			grg.setColor(Color.WHITE);
			grg.fillRect(0, 0, graphs.getWidth(), graphs.getHeight());


			//draw density (average front and back)
			//normalize density for drawing
			min = 100;
			max = 0.0;
			for(Boat bn : b){
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
			for(Boat bn: b){
				dpix = posy +graphHeight - (int)(graphHeight*((density[bn.index]-min)/(max-min)));
				grg.fillOval((int)(bn.x/R * graphWidth), dpix-2, 5, 5);
			}
			grg.drawString("density(boats/cm)", 50, posy - 2);
			grg.drawString("max = "+String.format("%.2f",max), 400, posy - 2);
			grg.drawString("min = "+String.format("%.2f",min), 600, posy - 2);
			
			//draw change in surface tension
			posy += graphHeight+15;
			
			//normalize dT for drawing
			min = 100;
			max = 0.0;
			for(Boat bn : b){
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
			for(Boat bn: b){
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
			min = 100;
			max = -100;
			for(Boat bn : b){
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
			for(Boat bn: b){
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
			min = 100;
			max = -1;
			for(Boat bn : b){
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
			for(Boat bn: b){
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
			min = 100000;
			max = -10;
			double totalfl = 0.0;
			for(Boat bn : b){
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
			for(Boat bn: b){
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
	}

	//////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	/**
     *  This is a panel for displaying and modifying parameters
     */
	class PPanel extends JPanel {
		// there are lots of labels and fields
        JLabel NLabel, RLabel, mLabel, visLabel, LLabel, lLabel, rcLabel, DLabel, kLabel, aLabel, bLabel;
        public JTextField NField, RField, mField, visField, LField, lField, rcField, DField, kField, aField, bField;
		double[] params;
		
		public PPanel() {
			super();
			params = model.getParameters();
			
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
			DLabel = new JLabel("diffusion constant(cm/s) = ");
			DLabel.setHorizontalAlignment(JLabel.TRAILING);
			kLabel = new JLabel("evaporation constant(1/s) = ");
			kLabel.setHorizontalAlignment(JLabel.TRAILING);
			aLabel = new JLabel("supply rate = ");
			aLabel.setHorizontalAlignment(JLabel.TRAILING);
			bLabel = new JLabel("beta = ");
			bLabel.setHorizontalAlignment(JLabel.TRAILING);
			
			// insert the initial parameters
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
				params[0] = Double.valueOf(RField.getText());
				params[1] = Double.valueOf(mField.getText());
				params[2] = Double.valueOf(visField.getText());
				params[3] = Double.valueOf(LField.getText());
				params[4] = Double.valueOf(lField.getText());
				params[5] = Double.valueOf(rcField.getText());
				params[6] = Double.valueOf(DField.getText());
				params[7] = Double.valueOf(kField.getText());
				params[8] = Double.valueOf(aField.getText());
				params[9] = Double.valueOf(bField.getText());
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