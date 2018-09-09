import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.util.*;
import java.io.*;
import javax.imageio.ImageIO;


/**
 * this is an external window for showing a plot
 */
class PlotWindow extends JPanel implements ActionListener{
	int width, height;
	JFrame frame;
	Plot2D plot;
	JButton saveButton;
	
	public PlotWindow(Plot2D theplot){
		super();
		plot = theplot;
		width = plot.width;
		height = plot.height;
		
		setPreferredSize(new Dimension(width, height));
		
		saveButton = new JButton("save plot");
		saveButton.addActionListener(this);
		
		frame = new JFrame("a plot");
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		
		Container pane = frame.getContentPane();
		pane.setLayout(new BorderLayout());
		pane.add(this, BorderLayout.CENTER);
		pane.add(saveButton, BorderLayout.SOUTH);
		
		pane.validate();
		//Display the window.
		frame.pack();
		frame.setVisible(true);
	}
	
	// for the actionlistener
	public void actionPerformed(ActionEvent e){
		if(e.getSource() == saveButton){
			writeImage();
		}  
	}
	
	//save an image in a file
	public void writeImage(){
		if(plot != null){
			File imageFile = new File("clusterImages");
			if(!imageFile.exists()){
				imageFile.mkdir();
			}
			imageFile = new File("clusterImages/test"+String.valueOf(System.currentTimeMillis())+".bmp");	//for everything else
			try{
				ImageIO.write(plot, "bmp", imageFile);
				JOptionPane.showMessageDialog(frame, "successfully recorded image\nFile: "+imageFile.getName());
			}catch(Exception e){
				System.err.println("Error: " + e.getMessage());
				JOptionPane.showMessageDialog(frame, "file writing error");
			}
		}
	}
	
	// ye olde painting method
	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		
		g.drawImage(plot, 0, 0, this);
		
	}
	
}
