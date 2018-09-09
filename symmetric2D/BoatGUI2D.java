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
class BoatGUI2D implements ActionListener{
    // These are the parts of the GUI
	BoatModel2D model;

	BoatsWorker worker;

	JFrame frame;
	BPanel boatsPanel;
	PPanel paramPanel;
	JPanel cntrlPanel, rightPanel, bottomPanel;
	
	BufferedImage experiment, graphs, trajectories, cview, stview;
	BufferedImage spaceTime;
	int imageType;
	String[] imageNames;
	boolean trajectoriesSet;
	
	int simWidth, simHeight;
	double guiTimeStart, guiTimeEnd;

	JMenuBar menuBar;
	JMenu modelMenu, windowMenu, simMenu;
    JMenuItem infoItem, fRateItem, fastItem, recordImageItem, sim1Item, sim2Item, sim3Item, sim4Item, sim5Item, sim6Item;
	ButtonGroup windowGroup;
	JRadioButtonMenuItem bigWindowItem, smallWindowItem;

    JButton runButton, timeButton, adButton, rmButton, homoButton, randButton, statsButton, discButton;
	JButton paramButton, resetButton, stopButton, displayButton, kickButton, writeButton;

	JTextField backgField;

	boolean showc;
	long compytime;
	int frameRate;
	boolean computeFast, recordData, recordEverything;

	File sim1File, sim2File, imageFile, everythingFile;
	FileWriter writer, nondimWriter, omniWriter, nondimomniWriter;
	FileOutputStream outStream;
	DataOutputStream dataOut;

	/**
	 * The constructor initializes all the model parts and prepares
	 * everything to run.
	 */
	public BoatGUI2D(){
		frameRate = 100;
		showc = false;
		computeFast = false;
		recordData = false;
		recordEverything = false;
		
		model = new BoatModel2D();

		worker = new BoatsWorker();

		simWidth = 800;
		simHeight = 600;
		imageType = 0;
		imageNames = new String[]{"experiment", "trajectories", "graphs", "c", "s. tension"};
		
		spaceTime = new BufferedImage(simWidth, simHeight, BufferedImage.TYPE_INT_ARGB);
		experiment = new BufferedImage(simWidth, simHeight, BufferedImage.TYPE_INT_ARGB);
		graphs = new BufferedImage(simWidth, simHeight, BufferedImage.TYPE_INT_ARGB);
		trajectories = new BufferedImage(simWidth, simHeight, BufferedImage.TYPE_INT_ARGB);
		cview = new BufferedImage(simWidth, simHeight, BufferedImage.TYPE_INT_ARGB);
		stview = new BufferedImage(simWidth, simHeight, BufferedImage.TYPE_INT_ARGB);
		trajectoriesSet = false;

		guiTimeStart = 0.0;
		guiTimeEnd = model.getInterval();

		frame = new JFrame("2D boats simulation");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		//Create the menu bar.
        menuBar = new JMenuBar();
        menuBar.setOpaque(true);
		//create some menus
		modelMenu = new JMenu("model");
		windowMenu = new JMenu("window");
		simMenu = new JMenu("simulations");
		//create some menu items
        infoItem = new JMenuItem("model info");
		infoItem.addActionListener(this);

		sim1Item = new JMenuItem("compute stats");
		sim1Item.addActionListener(this);
		sim2Item = new JMenuItem("vary N comp. stats");
		sim2Item.addActionListener(this);
		sim3Item = new JMenuItem("vary vis comp. v data");
		sim3Item.addActionListener(this);
		sim4Item = new JMenuItem("(not ready)");
		sim4Item.addActionListener(this);
		sim5Item = new JMenuItem("(not ready)");
		sim5Item.addActionListener(this);
		sim6Item = new JMenuItem("(not ready)");
		sim6Item.addActionListener(this);

		fRateItem = new JMenuItem("frame rate");
		fRateItem.addActionListener(this);
		fastItem = new JMenuItem("compute fast");
		fastItem.addActionListener(this);
		recordImageItem = new JMenuItem("save image");
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
		
		menuBar.add(modelMenu);
		menuBar.add(windowMenu);
		menuBar.add(simMenu);

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
		homoButton = new JButton("one step");
		homoButton.addActionListener(this);
		randButton = new JButton("randomize");
		randButton.addActionListener(this);
		statsButton = new JButton("stats");
		statsButton.addActionListener(this);
		discButton = new JButton("set dx, dt");
		discButton.addActionListener(this);
		paramButton = new JButton("set parameters");
		paramButton.addActionListener(this);
		resetButton = new JButton("reset");
		resetButton.addActionListener(this);
		stopButton = new JButton("stop");
		stopButton.addActionListener(this);
		displayButton = new JButton("change image");
		displayButton.addActionListener(this);
		kickButton = new JButton("kick a boat");
		kickButton.addActionListener(this);
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
		cntrlPanel.add(stopButton);
		cntrlPanel.add(timeButton);
		cntrlPanel.add(discButton);
		cntrlPanel.add(adButton);
        cntrlPanel.add(rmButton);
		cntrlPanel.add(randButton);
		cntrlPanel.add(statsButton);
		cntrlPanel.add(writeButton);
		cntrlPanel.add(new JLabel(" "));
		cntrlPanel.add(displayButton);
		cntrlPanel.add(homoButton);
		cntrlPanel.add(resetButton);
		cntrlPanel.add(paramButton);

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
			if(model.N > 0 && (worker.isCancelled() || !worker.getState().equals(SwingWorker.StateValue.STARTED))){
				trajectoriesSet = false;
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
			model.addBoat();
			paramPanel.NField.setText(String.valueOf(model.N));
			paramPanel.repaint();
			boatsPanel.repaint();
		} else if (e.getSource() == rmButton){
			if(worker != null){
				worker.cancel(true);
			}
			if(model.N>0){
				model.removeBoat();
			}
			paramPanel.NField.setText(String.valueOf(model.N));
			paramPanel.repaint();
			boatsPanel.repaint();
		} else if (e.getSource() == homoButton){
			model.homogenize();
			if(model.N > 0 && !worker.getState().equals(SwingWorker.StateValue.STARTED)){
				model.timeStep();
			}
			boatsPanel.repaint();
		} else if (e.getSource() == randButton){
			model.randomize();
			boatsPanel.repaint();
		} else if (e.getSource() == statsButton){
			if(!model.computeStats){
				if(model.N > 0 && (worker.isCancelled() || !worker.getState().equals(SwingWorker.StateValue.STARTED))){
					model.startRecording();
				}
			}else{
				model.stopRecording();
			}
			boatsPanel.repaint();
		} else if (e.getSource() == discButton){
			setDiscretization();
			boatsPanel.repaint();
		} else if (e.getSource() == paramButton){
			paramPanel.newParameters();
			boatsPanel.repaint();
		} else if (e.getSource() == resetButton){
			reset();
			paramPanel.NField.setText(String.valueOf(model.N));
			paramPanel.repaint();
		} else if (e.getSource() == stopButton){
			if(worker != null){
				worker.cancel(true);
			}
		} else if (e.getSource() == writeButton){
			if(worker != null){
				worker.cancel(true);
			}
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
		} else if (e.getSource() == kickButton){
			model.kickABoat(0,1000.0);
		} else if (e.getSource() == infoItem){
			JOptionPane.showMessageDialog(frame, "Hold on a moment, I will explain.");
		} else if (e.getSource() == fRateItem){
			frameRate = Integer.valueOf((String)JOptionPane.showInputDialog(frame, "Input an integer(~100)\n (this is the frame rate in model time, NOT real time)"));
		} else if (e.getSource() == fastItem){
			computeFast = !computeFast;
		} else if (e.getSource() == recordImageItem){
			if(!worker.getState().equals(SwingWorker.StateValue.STARTED)){
				double[] pars = model. getParameters();
				writeImage(boatsPanel.image, "images/N"+String.valueOf(model.N)+"vis"+String.valueOf((int)(10000*pars[3]))+"K"+String.valueOf((int)pars[8])+"T"+String.valueOf((int)model.getTime())+".bmp");
			}
			
		} else if (e.getSource() == bigWindowItem){
			boatsPanel.setSize(900, 800);
			simHeight = 600;
			reset();
			guiTimeStart = model.getTime();
			guiTimeEnd = guiTimeStart+model.getInterval();
			frame.pack();
			frame.repaint();
		} else if (e.getSource() == smallWindowItem){
			boatsPanel.setSize(900, 600);
			simHeight = 400;
			reset();
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
			runSim();
		} else if (e.getSource() == sim5Item){
			runSim();
		} else if (e.getSource() == sim6Item){
			runSim();
		}
	}
	
	/////////////////////////////////////////////////////////////////////////////////////////
	// beginning of simulations //
	/////////////////////////////////////////////////////////////////////////////////////////

	/**
	 * This computes the statistics for one run and records them
	 *
	 */
	public void runSim1(){
		if(JOptionPane.showConfirmDialog(frame, "Are all the parameters set?\nThis will begin running now.", "choose one", JOptionPane.YES_NO_OPTION) == JOptionPane.NO_OPTION){
			return;
		}

		if(model.N > 0 && (worker.isCancelled() || !worker.getState().equals(SwingWorker.StateValue.STARTED))){
			// keep going
		}else{
			JOptionPane.showMessageDialog(frame, "Wait, I wasn't ready. Try again.");
			return;
		}

		//write it to a file in columns
		sim1File = new File("data");
		if(!sim1File.exists()){
			sim1File.mkdir();
		}
		sim1File = new File("data/sim1stuff");
		if(!sim1File.exists()){
			sim1File.mkdir();
		}
		sim1File = new File("data/sim1stuff/statsN"+String.valueOf(model.N)+"R"+String.valueOf((int)model.R)+"vis"+String.valueOf((int)(1000*model.vis))+"t"+String.valueOf(System.currentTimeMillis()));
		File nondimFile = new File("data/sim1stuff/NONDIMstatsN"+String.valueOf(model.N)+"R"+String.valueOf((int)model.R)+"vis"+String.valueOf((int)(1000*model.vis))+"t"+String.valueOf(System.currentTimeMillis()));
		try{
			writer = new FileWriter(sim1File);
			writer.write("# time step, mean disp, mean sqr disp, ave min dist, ave v\n");
			nondimWriter = new FileWriter(nondimFile);
			nondimWriter.write("# time step, mean disp, mean sqr disp, ave min dist, ave v\n");
		}catch(Exception e){
			System.err.println("Error: " + e.getMessage());
			JOptionPane.showMessageDialog(frame, "file writing error");
		}
		
		// this is for the everything file
		everythingFile = new File("data/everything");
		if(!everythingFile.exists()){
			everythingFile.mkdir();
		}
		everythingFile = new File("data/everything/alldataN"+String.valueOf(model.N)+"t"+String.valueOf(System.currentTimeMillis()));
		try{
			outStream = new FileOutputStream(everythingFile);
			dataOut = new DataOutputStream(outStream);
		}catch(Exception e){
			System.err.println("Error: " + e.getMessage());
			JOptionPane.showMessageDialog(frame, "error creating everythingfile");
		}
		recordEverything = true;

		// now run it once whilst recording data
		trajectoriesSet = false;
		model.startRecording();
		computeFast = true;
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
		model.stopRecording();
		computeFast = false;
		
		recordEverything = false;
		try{
			outStream.close();
		}catch(Exception e){
			System.err.println("Error: " + e.getMessage());
			JOptionPane.showMessageDialog(frame, "file writing error");
		}

		// now write the data to the file
		try{
			for(int i=0; i<model.meanDisp.length; i++){
				writer.write(String.valueOf(i*4*model.dt)+" "+String.valueOf(model.meanDisp[i])+" "+String.valueOf(model.meansqDisp[i])+" "+String.valueOf(model.aveminDist[i])+" "+String.valueOf(model.avev[i])+"\n");
				nondimWriter.write(String.valueOf(i*4*model.dt*model.D/(model.L*model.L))+" "+String.valueOf(model.meanDisp[i]/model.L)+" "+String.valueOf(model.meansqDisp[i]/(model.L*model.L))+" "+String.valueOf(model.aveminDist[i]/model.L)+" "+String.valueOf(model.avev[i]*model.L/model.D)+"\n");
			}

		}catch(Exception e){
			System.err.println("Error: " + e.getMessage());
			JOptionPane.showMessageDialog(frame, "file writing error");
		}
		
		try{
			writer.close();
			nondimWriter.close();
			JOptionPane.showMessageDialog(frame, "successfully recorded data");
		}catch(Exception e){
			System.err.println("Error: " + e.getMessage());
			JOptionPane.showMessageDialog(frame, "file writing error");
		}
	}
	
	/**
	 * This computes the statistics for several values of N and records them
	 *
	 */
	public void runSim2(){
		int minN = model.N;
		int intervalN, stepsN;
		double relaxT, recordT;
		long totalStartTime;
		
		if(model.N > 0 && (worker.isCancelled() || !worker.getState().equals(SwingWorker.StateValue.STARTED))){
			// keep going
		}else{
			JOptionPane.showMessageDialog(frame, "Wait, I wasn't ready. Try again.");
			return;
		}
		
		if(JOptionPane.showConfirmDialog(frame, "Are all the parameters set?\nIs N at the minimum value?", "choose one", JOptionPane.YES_NO_OPTION) == JOptionPane.NO_OPTION){
			return;
		}
		try{
			intervalN = Integer.valueOf((String)JOptionPane.showInputDialog(frame, "Input N step size"));
			stepsN = Integer.valueOf((String)JOptionPane.showInputDialog(frame, "Input number of steps"));
			relaxT = Double.valueOf((String)JOptionPane.showInputDialog(frame, "Input time interval to run each simulation before recording"));
			recordT = Double.valueOf((String)JOptionPane.showInputDialog(frame, "Input time interval to record"));
		}catch(Exception e){
			return;
		}
		
		computeFast = true;
		
		// this is for the omnifile
		File omniFile = new File("data");
		if(!omniFile.exists()){
			omniFile.mkdir();
		}
		omniFile = new File("data/sim2stuff");
		if(!omniFile.exists()){
			omniFile.mkdir();
		}
		omniFile = new File("data/sim2stuff/statsomniN"+String.valueOf(model.N)+"R"+String.valueOf((int)model.R)+"vis"+String.valueOf((int)(1000*model.vis))+"t"+String.valueOf(System.currentTimeMillis()));
		File nondimomniFile = new File("data/sim2stuff/NONDIMstatsomniN"+String.valueOf(model.N)+"R"+String.valueOf((int)model.R)+"vis"+String.valueOf((int)(1000*model.vis))+"t"+String.valueOf(System.currentTimeMillis()));
		try{
			omniWriter = new FileWriter(omniFile);
			omniWriter.write("# N, time step, mean disp, mean sqr disp, ave min dist, ave v\n");
			nondimomniWriter = new FileWriter(nondimomniFile);
			nondimomniWriter.write("# N, time step, mean disp, mean sqr disp, ave min dist, ave v\n");
		}catch(Exception e){
			System.err.println("Error: " + e.getMessage());
			JOptionPane.showMessageDialog(frame, "error creating file");
		}
		
		// for each value of N: run initial period, run recording period, write data to separate file plus one omnifile
		totalStartTime = System.currentTimeMillis();
		for(int step=0; step<stepsN; step++){
			// reset and set up parameters
			reset();
			model.stopRecording();
			recordEverything = false;
			
			paramPanel.NField.setText(String.valueOf(minN+step*intervalN));
			paramButton.doClick();
			
			// run for initial period
			model.setInterval(relaxT);
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
			
			// run for recording period
			model.setInterval(recordT);
			model.startRecording();
			
			// this is for the everything file
			everythingFile = new File("data/everything");
			if(!everythingFile.exists()){
				everythingFile.mkdir();
			}
			everythingFile = new File("data/everything/alldataN"+String.valueOf(model.N)+"t"+String.valueOf(System.currentTimeMillis()));
			try{
				outStream = new FileOutputStream(everythingFile);
				dataOut = new DataOutputStream(outStream);
			}catch(Exception e){
				System.err.println("Error: " + e.getMessage());
				JOptionPane.showMessageDialog(frame, "error creating everythingfile");
			}
			recordEverything = true;
			
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
			
			recordEverything = false;
			
			// make a new file
			sim1File = new File("data/sim2stuff/statsN"+String.valueOf(model.N)+"R"+String.valueOf((int)model.R)+"vis"+String.valueOf((int)(1000*model.vis))+"t"+String.valueOf(System.currentTimeMillis()));
			File nondimFile = new File("data/sim2stuff/NONDIMstatsN"+String.valueOf(model.N)+"R"+String.valueOf((int)model.R)+"vis"+String.valueOf((int)(1000*model.vis))+"t"+String.valueOf(System.currentTimeMillis()));
			try{
				writer = new FileWriter(sim1File);
				writer.write("# time step, mean disp, mean sqr disp, ave min dist, ave v\n");
				nondimWriter = new FileWriter(nondimFile);
				nondimWriter.write("# time step, mean disp, mean sqr disp, ave min dist, ave v\n");
			}catch(Exception e){
				System.err.println("Error: " + e.getMessage());
				JOptionPane.showMessageDialog(frame, "error creating file");
			}
			
			// now write the data to the file
			try{
				for(int i=0; i<model.meanDisp.length; i++){
					omniWriter.write(String.valueOf(model.N)+" "+String.valueOf(i*4*model.dt)+" "+String.valueOf(model.meanDisp[i])+" "+String.valueOf(model.meansqDisp[i])+" "+String.valueOf(model.aveminDist[i])+" "+String.valueOf(model.avev[i])+"\n");
					nondimomniWriter.write(String.valueOf(model.N)+" "+String.valueOf(i*4*model.dt*model.D/(model.L*model.L))+" "+String.valueOf(model.meanDisp[i]/model.L)+" "+String.valueOf(model.meansqDisp[i]/(model.L*model.L))+" "+String.valueOf(model.aveminDist[i]/model.L)+" "+String.valueOf(model.avev[i]*model.L/model.D)+"\n");
					
					writer.write(String.valueOf(i*4*model.dt)+" "+String.valueOf(model.meanDisp[i])+" "+String.valueOf(model.meansqDisp[i])+" "+String.valueOf(model.aveminDist[i])+" "+String.valueOf(model.avev[i])+"\n");
					nondimWriter.write(String.valueOf(i*4*model.dt*model.D/(model.L*model.L))+" "+String.valueOf(model.meanDisp[i]/model.L)+" "+String.valueOf(model.meansqDisp[i]/(model.L*model.L))+" "+String.valueOf(model.aveminDist[i]/model.L)+" "+String.valueOf(model.avev[i]*model.L/model.D)+"\n");
				}
				
			}catch(Exception e){
				System.err.println("Error: " + e.getMessage());
				JOptionPane.showMessageDialog(frame, "error writing data");
			}
			
			try{
				writer.close();
				nondimWriter.close();
			}catch(Exception e){
				System.err.println("Error: " + e.getMessage());
				JOptionPane.showMessageDialog(frame, "file writing error");
			}
			try{
				outStream.close();
			}catch(Exception e){
				System.err.println("Error: " + e.getMessage());
				JOptionPane.showMessageDialog(frame, "file writing error");
			}
			
			System.out.println("finished step "+String.valueOf(step+1)+" of "+String.valueOf(stepsN));
		}
		
		try{
			omniWriter.close();
			nondimomniWriter.close();
		}catch(Exception e){
			System.err.println("Error: " + e.getMessage());
			JOptionPane.showMessageDialog(frame, "file writing error");
		}
		
		model.stopRecording();
		computeFast = false;
		System.out.println("simulation completed in "+String.valueOf((double)((System.currentTimeMillis()-totalStartTime)/6000)/10)+" minutes");
		
	}
	
	/**
	 * This computes the statistics for varying viscosity (ONE RUN) and records them
	 *
	 */
	public void runSim3(){
		double minvis = model.vis;
		double interval;
		int steps;
		double relaxT, recordT;
		long totalStartTime;
		
		if(model.N > 0 && (worker.isCancelled() || !worker.getState().equals(SwingWorker.StateValue.STARTED))){
			// keep going
		}else{
			JOptionPane.showMessageDialog(frame, "Wait, I wasn't ready. Try again.");
			return;
		}
		
		if(JOptionPane.showConfirmDialog(frame, "Are all the parameters set?\nAre they at the minimum value?", "choose one", JOptionPane.YES_NO_OPTION) == JOptionPane.NO_OPTION){
			return;
		}
		try{
			interval = Double.valueOf((String)JOptionPane.showInputDialog(frame, "Input viscosity step size"));
			steps = Integer.valueOf((String)JOptionPane.showInputDialog(frame, "Input number of steps"));
			relaxT = Double.valueOf((String)JOptionPane.showInputDialog(frame, "Input initial relaxation period"));
			recordT = Double.valueOf((String)JOptionPane.showInputDialog(frame, "Input time interval between measurements"));
		}catch(Exception e){
			return;
		}
		
		computeFast = true;
		
		// make a new file
		sim1File = new File("data");
		if(!sim1File.exists()){
			sim1File.mkdir();
		}
		sim1File = new File("data/sim3stuff");
		if(!sim1File.exists()){
			sim1File.mkdir();
		}
		sim1File = new File("data/sim3stuff/statsN"+String.valueOf(model.N)+"R"+String.valueOf((int)model.R)+"vis"+String.valueOf((int)(1000*model.vis))+"t"+String.valueOf(System.currentTimeMillis()));
		File nondimFile = new File("data/sim3stuff/NONDIMstatsN"+String.valueOf(model.N)+"R"+String.valueOf((int)model.R)+"vis"+String.valueOf((int)(1000*model.vis))+"t"+String.valueOf(System.currentTimeMillis()));
		try{
			writer = new FileWriter(sim1File);
			writer.write("# vis, mean(v), meansq(v), var(v)\n");
			nondimWriter = new FileWriter(nondimFile);
			nondimWriter.write("# vis, mean(v), meansq(v), var(v)\n");
		}catch(Exception e){
			System.err.println("Error: " + e.getMessage());
			JOptionPane.showMessageDialog(frame, "error creating file");
		}
		
		// for each value of vis: run initial period, run recording period, compute data
		totalStartTime = System.currentTimeMillis();
		
		int tmpN = model.N;
		reset();
		paramPanel.NField.setText(String.valueOf(tmpN));
		paramButton.doClick();
		
		// run for initial period
		model.setInterval(relaxT);
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
		model.setInterval(recordT);
		
		for(int step=0; step<steps; step++){
			// set up parameters
			paramPanel.visField.setText(String.valueOf(minvis+step*interval));
			paramButton.doClick();
			
			// run for recording period
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
			
			// compute the data
			double meanvx = 0.0;
			double meanvy = 0.0;
			double meansqv = 0.0;
			double varv = 0.0;
			double tmpvx, tmpvy;
			for(int i=0; i<model.N; i++){
				tmpvx = model.b.get(i).vx;
				tmpvy = model.b.get(i).vy;
				
				meanvx += tmpvx;
				meanvy += tmpvy;
				meansqv += tmpvx*tmpvx + tmpvy*tmpvy;
			}
			meanvx = meanvx/model.N;
			meanvy = meanvy/model.N;
			meanvx = Math.sqrt(meanvx*meanvx+meanvy*meanvy);
			meansqv = Math.sqrt(meansqv/model.N);
			varv = meansqv*meansqv - meanvx*meanvx;
			
			// now write the data to the file
			try{
				writer.write(String.valueOf(model.vis)+" "+String.valueOf(meanvx)+" "+String.valueOf(meansqv)+" "+String.valueOf(varv)+"\n");
				nondimWriter.write(String.valueOf(model.vis*model.L*model.L/model.m/model.D)+" "+String.valueOf(meanvx*model.L/model.D)+" "+String.valueOf(meansqv*model.L/model.D)+" "+String.valueOf(varv*model.L/model.D*model.L/model.D)+"\n");
			}catch(Exception e){
				System.err.println("Error: " + e.getMessage());
				JOptionPane.showMessageDialog(frame, "error writing data");
			}
			
			System.out.println("finished step "+String.valueOf(step+1)+" of "+String.valueOf(steps));
		}
		
		try{
			writer.close();
			nondimWriter.close();
		}catch(Exception e){
			System.err.println("Error: " + e.getMessage());
			JOptionPane.showMessageDialog(frame, "file writing error");
		}
		
		computeFast = false;
		System.out.println("simulation completed in "+String.valueOf((double)((System.currentTimeMillis()-totalStartTime)/6000)/10)+" minutes");
		
	}

	/**
	 * This is just a skeleton for future simulations
	 * 
	 */
	public void runSim(){
		// does nothing yet
	}
	
	
	//////////////////////////////////////////////////////////////////////////////////////
	// end of simulations //
	//////////////////////////////////////////////////////////////////////////////////////

	public void computeData(){
		
		// nothing yet
	}
	
	// writes all recorded stats to a file in columns
	public void writeData(){
		// if no stats are recorded, write nothing?
		if(model.meanDisp.length < 2){
			JOptionPane.showMessageDialog(frame, "Nothing has been computed yet.");
			return;
		}
		
		// make the files and directories
		sim1File = new File("data");
		if(!sim1File.exists()){
			sim1File.mkdir();
		}
		sim1File = new File("data/stats");
		if(!sim1File.exists()){
			sim1File.mkdir();
		}
		sim1File = new File("data/stats/statsN"+String.valueOf(model.N)+"R"+String.valueOf((int)model.R)+"vis"+String.valueOf((int)(1000*model.vis))+"t"+String.valueOf(System.currentTimeMillis()));
		File nondimFile = new File("data/stats/NONDIMstatsN"+String.valueOf(model.N)+"R"+String.valueOf((int)model.R)+"vis"+String.valueOf((int)(1000*model.vis))+"t"+String.valueOf(System.currentTimeMillis()));
		try{
			writer = new FileWriter(sim1File);
			writer.write("# time step, mean disp, mean sqr disp, ave min dist, ave v\n");
			nondimWriter = new FileWriter(nondimFile);
			nondimWriter.write("# time step, mean disp, mean sqr disp, ave min dist, ave v\n");
		}catch(Exception e){
			System.err.println("Error: " + e.getMessage());
			JOptionPane.showMessageDialog(frame, "error creating file");
		}
		// now write the data to the file
		try{
			for(int i=0; i<model.meanDisp.length; i++){
				writer.write(String.valueOf(i*4*model.dt)+" "+String.valueOf(model.meanDisp[i])+" "+String.valueOf(model.meansqDisp[i])+" "+String.valueOf(model.aveminDist[i])+" "+String.valueOf(model.avev[i])+"\n");
				nondimWriter.write(String.valueOf(i*4*model.dt*model.D/(model.L*model.L))+" "+String.valueOf(model.meanDisp[i]/model.L)+" "+String.valueOf(model.meansqDisp[i]/(model.L*model.L))+" "+String.valueOf(model.aveminDist[i]/model.L)+" "+String.valueOf(model.avev[i]*model.L/model.D)+"\n");
			}
			
		}catch(Exception e){
			System.err.println("Error: " + e.getMessage());
			JOptionPane.showMessageDialog(frame, "error writing data");
		}
		// close it
		try{
			writer.close();
			nondimWriter.close();
			JOptionPane.showMessageDialog(frame, "successfully recorded data");
		}catch(Exception e){
			System.err.println("Error: " + e.getMessage());
			JOptionPane.showMessageDialog(frame, "file writing error");
		}
	}
	
	// writes a supplied image to a specified file
	public void writeImage(BufferedImage img, String name){
		if(img != null){
			imageFile = new File("images");
			if(!imageFile.exists()){
				imageFile.mkdir();
			}
			double[] pars = model. getParameters();
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
	
	public void setDiscretization(){
		int nxperrc, ntpersec;
		// first ask for the spatial discretization, then time
		try{
			nxperrc = Integer.valueOf((String)JOptionPane.showInputDialog(frame, "Input number of grid points over rc"));
			ntpersec = Integer.valueOf((String)JOptionPane.showInputDialog(frame, "Input number of time steps per second\n (suggested >"+String.valueOf((int)(model.D*nxperrc*nxperrc/(model.rc*model.rc)))+")"));
			model.setDisc(nxperrc, ntpersec);
		}catch(Exception e){
			System.err.println("Error: " + e.getMessage());
			return;
		}
		
	}
	
	public void reset(){
		worker.cancel(true);
		worker = new BoatsWorker();
		model.reset();
		
		guiTimeStart = model.getTime();
		guiTimeEnd = guiTimeStart+model.getInterval();
		
		spaceTime = new BufferedImage(simWidth, simHeight, BufferedImage.TYPE_INT_ARGB);
		experiment = new BufferedImage(simWidth, simHeight, BufferedImage.TYPE_INT_ARGB);
		graphs = new BufferedImage(simWidth, simHeight, BufferedImage.TYPE_INT_ARGB);
		trajectories = new BufferedImage(simWidth, simHeight, BufferedImage.TYPE_INT_ARGB);
		cview = new BufferedImage(simWidth, simHeight, BufferedImage.TYPE_INT_ARGB);
		stview = new BufferedImage(simWidth, simHeight, BufferedImage.TYPE_INT_ARGB);

		boatsPanel.repaint();
	}

	//////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////
	/**
	 * This is the worker thread which runs the model calculations
	 * in the background and publishes the results
	 */
	private class BoatsWorker extends SwingWorker<String, BufferedImage> {
		Integer dummy = new Integer(1);
		
		@Override
        protected String doInBackground() {
            double[] disc = model.getDisc();
			int totalFrames = 1;
			int stepsPerChunk = (int)disc[2];
			int extraSteps = 0;
			
			if(frameRate > 0){
				totalFrames = (int)(disc[0]*frameRate);
				stepsPerChunk = (int)(disc[2] / totalFrames);
				if(stepsPerChunk < 1){
					totalFrames = (int)disc[2];
					stepsPerChunk = 1;
					extraSteps = 1;
				}else{
					extraSteps = ((int)disc[2])%(stepsPerChunk*totalFrames);
					while(extraSteps >= stepsPerChunk){
						extraSteps -= stepsPerChunk;
						totalFrames++;
					}
				}
			}
			if(computeFast){
				totalFrames = 1;
				stepsPerChunk = (int)disc[2];
				extraSteps = 0;
			}
			try{
				if(recordEverything){
					double[] theParameters = model.getParameters();
					// Warning: if dataOut has not been initialized ,this will crash
					dataOut.writeInt(stepsPerChunk*totalFrames+extraSteps);
					dataOut.writeInt(model.N);
					dataOut.writeDouble(theParameters[0]);
					dataOut.writeDouble(theParameters[3]);
					dataOut.writeDouble(theParameters[4]);
					dataOut.writeDouble(theParameters[1]);
					dataOut.writeDouble(theParameters[2]);
					dataOut.writeDouble(theParameters[5]);
					dataOut.writeDouble(theParameters[6]);
					dataOut.writeDouble(theParameters[7]);
					dataOut.writeDouble(theParameters[8]);
					dataOut.writeDouble(theParameters[9]);
					dataOut.writeDouble(theParameters[10]);
					dataOut.writeDouble(theParameters[11]);
					dataOut.writeDouble(theParameters[12]);
					dataOut.writeDouble(theParameters[13]);
					dataOut.writeDouble(theParameters[14]);
					dataOut.writeDouble(theParameters[15]);
					// also, don't forget to close the writer
				}
			}catch(Exception e){
				System.err.println("Error1: " + e.getMessage());
				JOptionPane.showMessageDialog(frame, "error writing everything");
			}
			
			BufferedImage bgst = new BufferedImage(simWidth, simHeight, BufferedImage.TYPE_INT_ARGB);
			int[] xpix = new int[model.N];
			int tpix;
			double stTime = 0.0;
			ArrayList<Boat2D> bgb;
			
			// do the time steps
			for(int j=0; j<=totalFrames; j++){
				if(isCancelled()){
					System.out.println("got cancelled!");
					return "cancelled";
				}
				if(j == totalFrames){
					stepsPerChunk = extraSteps;
				}
				
				for(int i=0; i<stepsPerChunk; i++){
					model.timeStep();
					
					bgb = model.getBoats();
					tpix = (int)((bgst.getHeight()-1)*stTime/(guiTimeEnd-guiTimeStart)+1);
					xpix[0] = (int)((bgst.getWidth()-1)*bgb.get(0).x/disc[3])+1;
					bgst.setRGB(xpix[0], bgst.getHeight()-tpix, Color.BLUE.getRGB());
					for(int k=1; k<model.N; k++){
						xpix[k] = (int)((bgst.getWidth()-1)*bgb.get(k).x/disc[3])+1;
						bgst.setRGB(xpix[k], bgst.getHeight()-tpix, Color.RED.getRGB());
					}
					
					if(recordEverything){
						try{
							Boat2D bn;
							for(int ind=0; ind<model.N; ind++){
								bn = model.b.get(ind);
								dataOut.writeDouble(model.getTime());
								dataOut.writeDouble(bn.x);
								dataOut.writeDouble(bn.y);
								dataOut.writeDouble(bn.vx);
								dataOut.writeDouble(bn.vy);
							}
						}catch(Exception e){
							System.err.println("Error2: " + e.getMessage());
							JOptionPane.showMessageDialog(frame, "error writing everything");
						}
					}
					
					stTime += disc[1];
				}
				publish(bgst);
			}
			return "finished"; //just because they said to do it
        }
		
        @Override
        protected void process(java.util.List<BufferedImage> billy) {
			spaceTime.getRaster().setRect(billy.get(billy.size()-1).getRaster());
			boatsPanel.paintImmediately(0,0,boatsPanel.getWidth(),boatsPanel.getHeight());
        }
		
		@Override
		protected void done() {
			try{
				//get();
			}catch(Exception ignore){
				System.err.println("Error3: " + ignore.getMessage());
			}
			
			// this doesn't actually do anything now
			
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
		public BufferedImage image;

		public BPanel() {
			width = 900;
			height = 850;
			image = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
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
			Graphics ig = image.getGraphics();
			ig.setColor(Color.WHITE);
			ig.fillRect(0, 0, width, height);
			ig.setColor(Color.BLACK);

			// get some model dimensions
			double[] dims = model.getParameters();
			double R = dims[0];
			double L = dims[4];
			double rc = dims[5];

			int picposx = 20;
			int picposy = 20;

			//normalize c for drawing
			float minc = Float.MAX_VALUE-1;
			float maxc = -minc;
			for(int i=0; i<model.nx; i++){
				for(int j=0; j<model.ny; j++){
					if(model.c[i][j]>maxc)
						maxc = (float)model.c[i][j];
					if(model.c[i][j]<minc)
						minc = (float)model.c[i][j];
				}
			}
			if(minc>=maxc){
				minc = 0;
				maxc = 1;
			}
			float vminc = Float.MAX_VALUE-1;
			float vmaxc = -vminc;
			for(int i=0; i<model.nx; i++){
				for(int j=0; j<model.ny; j++){
					if(model.volc[i][j]>vmaxc)
						vmaxc = (float)model.volc[i][j];
					if(model.volc[i][j]<vminc)
						vminc = (float)model.volc[i][j];
				}
			}
			if(vminc>=vmaxc){
				vminc = 0;
				vmaxc = 1;
			}
			
			ig.drawString("camphor concentration, threshold = "+String.format("%.4f",1/dims[12]), 20, experiment.getHeight()+picposy+10);
			ig.drawString("max="+String.format("%.4f",maxc), 400, experiment.getHeight()+picposy+10);
			ig.drawString("min="+String.format("%.4f",minc), 550, experiment.getHeight()+picposy+10);

			//draw the big picture
			if(imageType == 0){
				drawExperiment();
				ig.drawImage(experiment, picposx, picposy, this);
				ig.setColor(Color.BLACK);
				ig.drawString("experiment view", simWidth/2+picposx-90, 15);
			}else if(imageType == 1){
				ig.drawImage(spaceTime, picposx, picposy, this);
				ig.setColor(Color.BLACK);
				ig.drawString("trajectories", simWidth/2+picposx-90, 15);
			}else if(imageType == 2){
				drawGraphs();
				ig.drawImage(graphs, picposx, picposy, this);
				ig.setColor(Color.BLACK);
				ig.drawString("some data", simWidth/2+picposx-90, 15);
			}else if(imageType == 3){
				drawCview(minc, maxc);
				ig.drawImage(cview, picposx, picposy, this);

				//draw the camphor concentration
				int cposy = picposy+cview.getHeight()+30;
				ig.setColor(Color.WHITE);
				ig.fillRect(picposx, cposy, cview.getWidth(), 100);
				ig.setColor(Color.BLACK);
				ig.drawRect(picposx, cposy, cview.getWidth(), 100);
				ig.setColor(Color.RED);
				int cpix = 0;
				for(int i=0; i<model.nx; i++){
					cpix = cposy + 100 - (int)(100*(model.c[i][model.ny/2]/maxc));
					ig.drawLine((int)(picposx+i*cview.getWidth()/model.nx), cpix, (int)(picposx+(i+1)*cview.getWidth()/model.nx), cpix);
				}
				ig.drawString(String.format("%.4f",maxc), simWidth+picposx+5, cposy+10);
				ig.drawString(String.format("%.4f",minc), simWidth+picposx+5, cposy+90);
				ig.drawString("surf.", simWidth+picposx+5, cposy+40);
				ig.setColor(new Color(0, 100, 0));
				int vcpix = 0;
				for(int i=0; i<model.nx; i++){
					vcpix = cposy + 100 - (int)(100*(model.volc[i][model.ny/2]/vmaxc));
					ig.drawLine((int)(picposx+i*cview.getWidth()/model.nx), vcpix, (int)(picposx+(i+1)*cview.getWidth()/model.nx), vcpix);
				}
				ig.drawString(String.format("%.4f",vmaxc), simWidth+picposx+5, cposy+20);
				ig.drawString(String.format("%.4f",vminc), simWidth+picposx+5, cposy+100);
				ig.drawString("vol.", simWidth+picposx+5, cposy+50);
				ig.setColor(Color.BLACK);
				ig.drawString("camphor concentration", simWidth/2+picposx-90, 15);
			}else if(imageType == 4){
				drawSTview(minc, maxc);
				ig.drawImage(stview, picposx, picposy, this);
				ig.setColor(Color.BLACK);
				ig.drawString("surface tension", simWidth/2+picposx-90, 15);
				
				// draw the color scale
				Color cola;
				double colval;
				double minst = model.maxdgamma/(model.beta*model.beta*maxc*maxc+1);
				double maxst = model.maxdgamma/(model.beta*model.beta*minc*minc+1);
				double dst = maxst-minst;
				int tmpwide = stview.getWidth();
				int cposy = picposy+stview.getHeight()+50;
				
				for(int i=0; i<tmpwide; i++){
					colval = (i*(model.maxdgamma/tmpwide)-minst)/dst;
					if(colval < 0){ colval = 0.0; }
					if(colval > 1){ colval = 1.0; }
					cola = Color.getHSBColor((float)(0.75-0.75*colval), (float)0.9, (float)1.0);
					ig.setColor(cola);
					ig.drawLine(picposx+i, cposy, picposx+i, cposy+20);
				}
				
				// draw the values
				ig.setColor(Color.BLACK);
				ig.drawRect(picposx, cposy, stview.getWidth(), 20);
				ig.drawString(String.format("%.3f",72-model.maxdgamma), picposx, cposy+35);
				ig.drawString("72.0", picposx+tmpwide-30, cposy+35);
				ig.drawString("surface tension color scale", picposx+tmpwide/2-80, cposy+35);
				int tmpx = picposx+(int)(tmpwide*minst/model.maxdgamma);
				ig.drawLine(tmpx, cposy-5, tmpx, cposy+20);
				ig.drawString(String.format("%.2f",72-model.maxdgamma+minst), tmpx+2, cposy-5);
				tmpx = picposx+(int)(tmpwide*maxst/model.maxdgamma);
				ig.drawLine(tmpx, cposy-5, tmpx, cposy+20);
				ig.drawString(String.format("%.2f",72-model.maxdgamma+maxst), tmpx+2, cposy-5);
			}

			//draw some useful stuff
			ig.setColor(Color.BLACK);
			
			ig.drawString("time", simWidth+picposx+5, (int)(experiment.getHeight()/2+picposy));
			ig.drawString("t = "+String.format("%.2f",model.getTime()), simWidth+picposx+5, (int)(experiment.getHeight()/2+picposy + 20));
			ig.drawString(String.format("%.2f",guiTimeEnd), simWidth+picposx+5, picposy+10);
			ig.drawString(String.format("%.2f",guiTimeStart), simWidth+picposx+5, experiment.getHeight()+picposy);
			
			g.drawImage(image, 0, 0, null);
		}
		
		public void drawExperiment(){
			double[] dims = model.getParameters();
			double R = dims[0];
			double W = dims[1];
			double L = dims[4];
			double rc = dims[5];
			Graphics2D exg = experiment.createGraphics();
			int exDia = experiment.getHeight()-40;
			double scale = (exDia)/(R/Math.PI + L);
			int boatDia = (int)(L*scale);
			int routew = (int)(W*scale);
			double angpercm = 2*Math.PI/model.R;
			double innerrad = cview.getHeight()/2 - 20 - model.W*scale;
			double outerrad = cview.getHeight()/2 - 20;
			
			exg.setColor(Color.WHITE);
			exg.fillRect(0, 0, experiment.getWidth(), experiment.getHeight());
			exg.setColor(Color.BLACK);
			exg.drawOval(20, 20, (int)(2*outerrad), (int)(2*outerrad));
			exg.drawOval(20+(int)(model.W*scale), 20+(int)(model.W*scale), (int)(2*innerrad), (int)(2*innerrad));
			exg.drawLine(50+exDia, 20, 50+exDia, 20+exDia);
			exg.drawLine(40+exDia, 20, 60+exDia, 20);
			exg.drawLine(40+exDia, 20+exDia, 60+exDia, 20+exDia);
			exg.drawString(String.format("%.2f", R/Math.PI)+"cm", 55+exDia, 10+exDia/2);
			for(Boat2D bn: model.b){
				int bxpos = (int) (experiment.getHeight()/2 + (innerrad + bn.y*scale)*Math.cos(bn.x*angpercm) - boatDia/2.0);
				int bypos = (int) (experiment.getHeight()/2 + (innerrad + bn.y*scale)*Math.sin(bn.x*angpercm) + boatDia/2.0);
				exg.fillOval(bxpos, experiment.getHeight()-bypos , boatDia, boatDia);
			}
			
		}
		
		public void drawTrajectories(){
			double[] dims = model.getParameters();
			double R = dims[0];
			Graphics2D tg = trajectories.createGraphics();
			
			int exDia = experiment.getHeight()-40;
			double scale = exDia/R;
			// draw the unchanging stuff
			if(!trajectoriesSet){
				tg.setColor(Color.WHITE);
				tg.fillRect(0, 0, trajectories.getWidth(), trajectories.getHeight());
				tg.setColor(Color.BLACK);
				tg.drawRect(20, 20, exDia, exDia);
				tg.drawLine(50+exDia, 20, 50+exDia, 20+exDia);
				tg.drawLine(40+exDia, 20, 60+exDia, 20);
				tg.drawLine(40+exDia, 20+exDia, 60+exDia, 20+exDia);
				tg.drawString(String.format("%.2f", R), 55+exDia, 10+exDia/2);
				trajectoriesSet = true;
			}
			
			int tDia = 3;
			int bxpos, bypos;
			tg.setColor(Color.BLACK);
			Boat2D bn;
			for(int i=0; i<model.N; i++){
				bn = model.b.get(i);
				bxpos = (int)(20 + (bn.x)*scale);
				bypos = (int)(20 + (bn.y)*scale);
				tg.fillOval(bxpos, trajectories.getHeight() - bypos , tDia, tDia);
			}
		}
		
		public void drawGraphs(){
			Graphics2D grg = graphs.createGraphics();
			grg.setColor(Color.WHITE);
			grg.fillRect(0,0,graphs.getWidth(),graphs.getHeight());
			if(model.computeStats){
				// draw 4 graphs showing meanDisp, meansqDisp, aveminDist, avev
				int xpt, ypt;
				double maxval;
				int wide = graphs.getWidth()-50;
				int high = graphs.getHeight()/4 - 4;

				//draw the axes
				grg.setColor(Color.BLACK);
				grg.drawLine(5,2,5,high+2);
				grg.drawLine(5,high+6,5,2*high+6);
				grg.drawLine(5,2*high+10,5,3*high+10);
				grg.drawLine(5,3*high+14,5,4*high+14);
				grg.drawLine(5,high+2,wide+5,high+2);
				grg.drawLine(5,2*high+6,wide+5,2*high+6);
				grg.drawLine(5,3*high+10,wide+5,3*high+10);
				grg.drawLine(5,4*high+14,wide+5,4*high+14);

				grg.setColor(Color.RED);
				// meanDisp
				maxval = 0.0;
				for(int i=0; i<model.meanDisp.length; i++){
					if(model.meanDisp[i] > maxval){
						maxval = model.meanDisp[i];
					}
				}
				for(int i=0; i<model.meanDisp.length; i++){
					xpt = (wide*i)/model.meanDisp.length+5;
					ypt = high+2 - (int)(model.meanDisp[i]/maxval*high);
					grg.drawRect(xpt,ypt,1,1);
				}
				grg.drawString(String.valueOf(maxval),5,12);
				grg.drawString("<X>",wide+5,12);

				// meansqDisp
				maxval = 0.0;
				for(int i=0; i<model.meansqDisp.length; i++){
					if(model.meansqDisp[i] > maxval){
						maxval = model.meansqDisp[i];
					}
				}
				for(int i=0; i<model.meansqDisp.length; i++){
					xpt = (wide*i)/model.meansqDisp.length+5;
					ypt = 2*high+6 - (int)(model.meansqDisp[i]/maxval*high);
					grg.drawRect(xpt,ypt,1,1);
				}
				grg.drawString(String.valueOf(maxval),5,high+16);
				grg.drawString("<x2>",wide+5,high+16);

				// aveminDist
				maxval = 0.0;
				for(int i=0; i<model.aveminDist.length; i++){
					if(model.aveminDist[i] > maxval){
						maxval = model.aveminDist[i];
					}
				}
				for(int i=0; i<model.aveminDist.length; i++){
					xpt = (wide*i)/model.aveminDist.length+5;
					ypt = 3*high+10 - (int)(model.aveminDist[i]/maxval*high);
					grg.drawRect(xpt,ypt,1,1);
				}
				grg.drawString(String.valueOf(maxval),5,2*high+20);
				grg.drawString("<rmin>",wide+5,2*high+20);

				// avev
				maxval = 0.0;
				for(int i=0; i<model.avev.length; i++){
					if(model.avev[i] > maxval){
						maxval = model.avev[i];
					}
				}
				for(int i=0; i<model.avev.length; i++){
					xpt = (wide*i)/model.avev.length+5;
					ypt = 4*high+14 - (int)(model.avev[i]/maxval*high);
					grg.drawRect(xpt,ypt,1,1);
				}
				grg.drawString(String.valueOf(maxval),5,3*high+24);
				grg.drawString("<|v|>",wide+5,3*high+24);

			}else{
				grg.setColor(Color.BLACK);
				grg.drawString("nothing here yet", 100, 100);
			}
		}
		
		public void drawCview(float minc, float maxc){
			Graphics2D cg = cview.createGraphics();
			cg.setColor(Color.WHITE);
			cg.fillRect(0,0,cview.getWidth(),cview.getHeight());

			Color color;
			//double scale = cview.getHeight()/model.R;
			//double scale = cview.getHeight()/(model.R/Math.PI);
			int exDia = experiment.getHeight()-40;
			double scale = (exDia)/(model.R/Math.PI + model.L/2.0);
			int boatDia = (int)(model.L*scale);
			int routew = (int)(model.W*scale);
			double angpercm = 2*Math.PI/model.R;
			double innerrad = cview.getHeight()/2 - 20 - model.W*scale;
			double outerrad = cview.getHeight()/2 - 20;
			
			cg.setColor(Color.BLACK);
			cg.drawOval(20, 20, (int)(2*outerrad), (int)(2*outerrad));
			cg.drawOval(20+(int)(model.W*scale), 20+(int)(model.W*scale), (int)(2*innerrad), (int)(2*innerrad));
			
			int pointwidth = cview.getWidth()/model.nx + 1;
			int xpix, ypix;
			for(int i=0; i<model.nx; i++){
				for(int j=0; j<model.ny; j++){
					//xpix = (int)(i*model.dx*scale);
					//ypix = (int)(j*model.dx*scale);
					///////////
					// paint it on a circle
					/////////////
					xpix = (int)(cview.getHeight()/2 + (innerrad + j*model.dx*scale)*Math.cos(i*model.dx*angpercm));
					ypix = (int)(cview.getHeight()/2 + (innerrad + j*model.dx*scale)*Math.sin(i*model.dx*angpercm));
					/////////////
					color = Color.getHSBColor((float)(0.75-0.75*((model.c[i][j]-minc)/(maxc-minc))), (float)1.0, (float)1.0);
					cg.setColor(color);
					cg.fillRect(xpix, cview.getHeight()-ypix, pointwidth, pointwidth);
				}
			}
			/*
			int boatDia = (int)(model.L*scale);
			int rDia = (int)(model.rc*2*scale);
			int bxpos, bypos, rxpos, rypos;
			cg.setColor(Color.BLACK);
			Boat2D bn;
			for(int i=0; i<model.N; i++){
				bn = model.b.get(i);
				bxpos = (int)((bn.x - model.L/2)*scale);
				bypos = (int)((bn.y + model.L/2)*scale);
				rxpos = (int)((bn.x - model.rc)*scale);
				rypos = (int)((bn.y + model.rc)*scale);

				cg.drawOval(bxpos, cview.getHeight() - bypos , boatDia, boatDia);
				cg.drawOval(rxpos, cview.getHeight() - rypos , rDia, rDia);
			}
			*/
			cg.setColor(Color.BLACK);
			for(Boat2D bn: model.b){
				int bxpos = (int) (cview.getHeight()/2 + (innerrad + bn.y*scale)*Math.cos(bn.x*angpercm) - boatDia/2.0);
				int bypos = (int) (cview.getHeight()/2 + (innerrad + bn.y*scale)*Math.sin(bn.x*angpercm) + boatDia/2.0);
				cg.drawOval(bxpos, cview.getHeight()-bypos , boatDia, boatDia);
			}
		}
		
		public void drawSTview(float minc, float maxc){
			Graphics2D stg = stview.createGraphics();
			stg.setColor(Color.WHITE);
			stg.fillRect(0,0,stview.getWidth(),stview.getHeight());
			
			Color color;
			double scale = stview.getHeight()/model.R;
			float minst = (float)(1.0/(model.beta*model.beta*maxc*maxc+1));
			float maxst = (float)(1.0/(model.beta*model.beta*minc*minc+1));
			int pointwidth = cview.getWidth()/model.nx + 1;
			int xpix, ypix;
			double st;
			for(int i=0; i<model.nx; i++){
				for(int j=0; j<model.ny; j++){
					xpix = (int)(i*model.dx*scale);
					ypix = (int)(j*model.dx*scale);
					st = (1.0/(model.beta*model.beta*model.c[i][j]*model.c[i][j]+1) - minst)/(maxst-minst);
					if(st > 1.0){st=1.0;}
					if(st <0.0){st=0.0;}
					color = Color.getHSBColor((float)(0.75-0.75*st), (float)0.9, (float)1.0);
					stg.setColor(color);
					stg.fillRect(xpix, stview.getHeight()-ypix, pointwidth, pointwidth);
				}
			}
			int boatDia = (int)(model.L*scale);
			int rDia = (int)(model.rc*2*scale);
			int bxpos, bypos, rxpos, rypos;
			double maxforce = 0.0;
			double tmpforce = 0.0;
			for(int i=0; i<model.N; i++){
				tmpforce = model.netx[i]*model.netx[i] + model.nety[i]*model.nety[i];
				if(maxforce < tmpforce){
					maxforce = tmpforce;
				}
			}
			maxforce = Math.sqrt(maxforce);
			stg.setColor(Color.BLACK);
			Boat2D bn;
			for(int i=0; i<model.N; i++){
				bn = model.b.get(i);
				bxpos = (int)((bn.x - model.L/2)*scale);
				bypos = (int)((bn.y + model.L/2)*scale);
				rxpos = (int)((bn.x - model.rc)*scale);
				rypos = (int)((bn.y + model.rc)*scale);
				
				stg.drawOval(bxpos, stview.getHeight() - bypos , boatDia, boatDia);
				stg.drawOval(rxpos, stview.getHeight() - rypos , rDia, rDia);
				stg.drawRect((int)(bn.x*scale)-1, stview.getHeight() - (int)(bn.y*scale)-1, 2, 2);
				stg.drawLine((int)(bn.x*scale), stview.getHeight() - (int)(bn.y*scale), (int)(bn.x*scale + model.netx[i]*boatDia/2/maxforce), stview.getHeight() - (int)(bn.y*scale + model.nety[i]*boatDia/2/maxforce));
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
        JLabel NLabel, RLabel, WLabel, mLabel, visLabel, LLabel, rcLabel, DLabel, kLabel, aLabel, vDLabel, vkLabel, vaLabel, bLabel, csatLabel, dgLabel, restLabel;
        public JTextField NField, RField, WField, mField, visField, LField, rcField, DField, kField, aField, vDField, vkField, vaField, bField, csatField, dgField, restField;

		public PPanel() {
			// there are lots of labels and fields
			NLabel = new JLabel("N = ");
            NLabel.setHorizontalAlignment(JLabel.TRAILING);
            RLabel = new JLabel("route length(cm) = ");
            RLabel.setHorizontalAlignment(JLabel.TRAILING);
			WLabel = new JLabel("route width(cm) = ");
            WLabel.setHorizontalAlignment(JLabel.TRAILING);
            mLabel = new JLabel("mass(g) = ");
            mLabel.setHorizontalAlignment(JLabel.TRAILING);
            visLabel = new JLabel("viscosity constant(g/s) = ");
            visLabel.setHorizontalAlignment(JLabel.TRAILING);
            LLabel = new JLabel("diameter of boat(cm) = ");
            LLabel.setHorizontalAlignment(JLabel.TRAILING);
            rcLabel = new JLabel("radius of camphor supply(cm) = ");
            rcLabel.setHorizontalAlignment(JLabel.TRAILING);
            DLabel = new JLabel("surface diffusion(cm*cm/s) = ");
            DLabel.setHorizontalAlignment(JLabel.TRAILING);
			vDLabel = new JLabel("volume diffusion(cm*cm/s) = ");
            vDLabel.setHorizontalAlignment(JLabel.TRAILING);
            kLabel = new JLabel("sublimation(1/s) = ");
            kLabel.setHorizontalAlignment(JLabel.TRAILING);
			vkLabel = new JLabel("surf/vol exchange(1/s) = ");
            vkLabel.setHorizontalAlignment(JLabel.TRAILING);
            aLabel = new JLabel("surf. supply rate = ");
            aLabel.setHorizontalAlignment(JLabel.TRAILING);
			vaLabel = new JLabel("vol. supply rate = ");
            vaLabel.setHorizontalAlignment(JLabel.TRAILING);
            bLabel = new JLabel("1/threshold = ");
            bLabel.setHorizontalAlignment(JLabel.TRAILING);
			csatLabel = new JLabel("c saturation = ");
            csatLabel.setHorizontalAlignment(JLabel.TRAILING);
			dgLabel = new JLabel("max dgamma = ");
            dgLabel.setHorizontalAlignment(JLabel.TRAILING);
			restLabel = new JLabel("restitution = ");
            restLabel.setHorizontalAlignment(JLabel.TRAILING);

			// insert the initial parameters
			double[] params =  model.getParameters();
			NField = new JTextField(String.valueOf(model.N), 6);
            NField.setEditable(true);
            RField = new JTextField(String.valueOf(params[0]), 6);
            RField.setEditable(true);
			WField = new JTextField(String.valueOf(params[1]), 6);
            WField.setEditable(true);
            mField = new JTextField(String.valueOf(params[2]), 6);
            mField.setEditable(true);
            visField = new JTextField(String.valueOf(params[3]), 6);
            visField.setEditable(true);
            LField = new JTextField(String.valueOf(params[4]), 6);
            LField.setEditable(true);
            rcField = new JTextField(String.valueOf(params[5]), 6);
            rcField.setEditable(true);
            DField = new JTextField(String.valueOf(params[6]), 6);
            DField.setEditable(true);
			vDField = new JTextField(String.valueOf(params[7]), 6);
            vDField.setEditable(true);
            kField = new JTextField(String.valueOf(params[8]), 6);
            kField.setEditable(true);
			vkField = new JTextField(String.valueOf(params[9]), 6);
            vkField.setEditable(true);
            aField = new JTextField(String.valueOf(params[10]), 6);
            aField.setEditable(true);
			vaField = new JTextField(String.valueOf(params[11]), 6);
            vaField.setEditable(true);
            bField = new JTextField(String.valueOf(params[12]), 6);
            bField.setEditable(true);
			csatField = new JTextField(String.valueOf(params[13]), 6);
            csatField.setEditable(true);
			dgField = new JTextField(String.valueOf(params[14]), 6);
            dgField.setEditable(true);
			restField = new JTextField(String.valueOf(params[15]), 6);
            restField.setEditable(true);

            // now set up the panel with the above items
            this.setLayout(new GridLayout(0,2));
			this.add(NLabel);
			this.add(NField);
            this.add(RLabel);
            this.add(RField);
			this.add(WLabel);
            this.add(WField);
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
			this.add(restLabel);
            this.add(restField);
		}

		public void newParameters(){
			try{
				double[] params = new double[16];
				params[0] = Double.valueOf(RField.getText());
				params[1] = Double.valueOf(WField.getText());
				params[2] = Double.valueOf(mField.getText());
				params[3] = Double.valueOf(visField.getText());
				params[4] = Double.valueOf(LField.getText());
				params[5] = Double.valueOf(rcField.getText());
				params[6] = Double.valueOf(DField.getText());
				params[7] = Double.valueOf(vDField.getText());
				params[8] = Double.valueOf(kField.getText());
				params[9] = Double.valueOf(vkField.getText());
				params[10] = Double.valueOf(aField.getText());
				params[11] = Double.valueOf(vaField.getText());
				params[12] = Double.valueOf(bField.getText());
				params[13] = Double.valueOf(csatField.getText());
				params[14] = Double.valueOf(dgField.getText());
				params[15] = Double.valueOf(restField.getText());
				
				model.setParameters(params);
				if(model.getN() != (int)Integer.valueOf(NField.getText())){
					model.setN((int)Integer.valueOf(NField.getText()));
				}
			} catch(Exception expt){
				JOptionPane.showMessageDialog(this, "Couldn't read parameters. Try again.");
			}
		}
	}

}