// Run -> Run File rather than Run Project

package atomviz;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.applet.*;

/**AtomViz: An interactive molecular dynamics applet to simulate atoms and
 * molecules in two dimensions.
 *
 * Author: Professor Paul Millett, University of Arkansas, pmillett@uark.edu
 *
 * Note: many aspects of this code were inspired by the software package,
 * MDApplet, written by Prof. Dan Schroeder, Physics Department, Weber State
 * University (see: http://physics.weber.edu/schroeder/software/mdapplet.html) */

/**Copyright @ 2016 Paul Millett & the Board of Trustees of the University of
 * Arkansas
 *
 * The terms under which this software, AtomViz, and associated documentation is
 * provided are as the following: 1. The Software is provided "as is," without
 * warranty of any kind, express or implied, including but not limited to the
 * warranties of merchantability, fitness for a particular purpose and
 * noninfringement. 2. In no event shall the Authors or Copyright Holder be
 * liable for any claim, damages or other liability, whether in an action of
 * contract, tort or otherwise, arising from, out of or in connection with
 * AtomViz or the use or other dealings in AtomViz. 3. The Authors or Copyright
 * Holder grant, free of charge, to any users the right to modify, copy, and
 * redistribute AtomViz, both within the user’s organization and externally,
 * subject to the following restrictions: a. The users agree not to charge for
 * the code itself. b. In any product based on AtomViz, the users agree to
 * acknowledge Prof. Millett. This acknowledgment shall appear in the project
 * documentation. c. The users agree to obey all U.S. Government restrictions
 * governing redistribution or export of The Software. d. The users agree to
 * reproduce any copyright notice which appears on this software on any copy or
 * modification of such made available to others.
 *
 * August 2016*/

public class AtomVizBase extends Applet implements Runnable {
    
    // ------------------------------------------------------------------------------
    // Data structures for AtomViz:
    // ------------------------------------------------------------------------------
    int pixelsPerUnit = 14;     // number of pixels in one distance unit
    int Nmax = 625;             // maximum number of particles
    int N = 625;                // actual number of particles
    double dt = 0.01;           // time step size
    double t = 0.0;             // simulation time
    double boxWidth = 49.8;     // simulation box width
    double wallStiffness = 50;  // 'spring constant' for bouncing off walls
    double gravity = 0.05;      // gravitational acceleration in the y-dir.
    double forceCutoff = 3;     // distance beyond which we don't bother to compute the force
    double forceCutoff2 = forceCutoff * forceCutoff;
    double initialKineticEnergy = 0.95;  // initial total kinetic energy
    double desiredKineticEnergy = 0.0;   // the specified current kinetic energy
    double maxKineticEnergy = 1.0;       // max total kinetic energy

    // here are the arrays of molecule positions, velocities, and accelerations:
    double[] x = new double[Nmax];
    double[] y = new double[Nmax];
    double[] vx = new double[Nmax];
    double[] vy = new double[Nmax];
    double[] ax = new double[Nmax];
    double[] ay = new double[Nmax];
    Color[] mColor = new Color[Nmax];    // colors for "random" setting

    // flags to show whether one is running base or their own class:
    boolean firstPrintOff = true;
    boolean derivedClassPos = true;
    boolean derivedClassVel = true;
    boolean derivedClassSep = true;

    // variables associated with drawing:
    int canvasWidth = 700;
    boolean running = false;        // true when simulation is running
    boolean reinit = false;         // true if "Restart" button has been pressed
    MDCanvas theCanvas;
    Button startButton;
    Button restartButton;
    Button addheatButton;
    Button cooldownButton;
    JProgressBar thermostat;
    JCheckBox gravitySwitch;
    JCheckBox slowmotionSwitch;

    // carefully constructed list of colors for indicating speeds:
    final Color[] speedColorList = {new Color(20, 0, 180), new Color(60, 0, 170), new Color(80, 0, 160), new Color(100, 0, 150),
        new Color(120, 0, 120), new Color(140, 0, 80), new Color(160, 0, 0), new Color(180, 0, 0),
        new Color(200, 0, 0), new Color(230, 0, 0), Color.RED, new Color(255, 60, 0),
        new Color(255, 90, 0), new Color(255, 120, 0), new Color(255, 150, 0), new Color(255, 180, 0),
        new Color(255, 210, 0), new Color(255, 230, 0), Color.YELLOW, new Color(255, 255, 120)};

    // ------------------------------------------------------------------------------
    // Exceedingly long init method, mostly to set up the GUI.  
    // ------------------------------------------------------------------------------
    @Override
    public void init() {
        // main elements and layout:
        setLayout(new BorderLayout());    // we'll use only CENTER, EAST, and SOUTH positions
        Color marginColor = Color.lightGray;
        Panel leftPanel = new Panel();
        leftPanel.setBackground(marginColor);
        add(leftPanel, BorderLayout.CENTER);
        theCanvas = new MDCanvas();
        leftPanel.add(theCanvas);
        Panel controlPanel = new Panel();
        add(controlPanel, BorderLayout.SOUTH);
        controlPanel.setBackground(marginColor);
        Panel ThermostatPanel = new Panel();
        leftPanel.add(ThermostatPanel);

        // here come the buttons...
        Panel startPanel = new Panel();
        startPanel.setLayout(new GridLayout(1, 2));
        controlPanel.add(startPanel);

        startButton = new Button("Start");
        startPanel.add(startButton);
        startButton.addActionListener((ActionEvent e) -> {
            running = !running;
            if (running) {
                startButton.setLabel("Pause");
            } else {
                startButton.setLabel("Resume");
            }
        });

        restartButton = new Button("Restart");
        startPanel.add(restartButton);
        restartButton.addActionListener((ActionEvent e) -> {
            reinit = true;
        });

        // thermostat bar...
        thermostat = new JProgressBar(0, 100);
        thermostat.setOrientation(1);
        Dimension thermSize = thermostat.getPreferredSize();
        thermSize.height = 700;
        thermostat.setPreferredSize(thermSize);
        ThermostatPanel.add(thermostat);
        thermostat.setValue((int) (initialKineticEnergy / maxKineticEnergy * 100));
        thermostat.setStringPainted(true);

        // temperature panel...
        Panel tempPanel = new Panel();
        tempPanel.setLayout(new GridLayout(1, 2));
        controlPanel.add(tempPanel);

        addheatButton = new Button("Add Heat");
        tempPanel.add(addheatButton);
        addheatButton.addActionListener((ActionEvent e) -> {
            desiredKineticEnergy *= 1.1;
            if (desiredKineticEnergy > maxKineticEnergy) {
                desiredKineticEnergy = maxKineticEnergy;
            }
        });

        cooldownButton = new Button("Cool Down");
        tempPanel.add(cooldownButton);
        cooldownButton.addActionListener((ActionEvent e) -> {
            desiredKineticEnergy *= 0.9;
        });

        // gravity checkbox:
        gravitySwitch = new JCheckBox("gravity");
        tempPanel.add(gravitySwitch);
        gravitySwitch.setSelected(true);

        // slow-motion checkbox:
        slowmotionSwitch = new JCheckBox("slow motion?");
        tempPanel.add(slowmotionSwitch);
        slowmotionSwitch.setSelected(false);

        // place the initial molecules:
        initializeParticles();

        // start a new thread to run the simulation:
        Thread thread = new Thread(this);
        thread.start();

    }    // end of init method

    // ------------------------------------------------------------------------------
    // Initialize particles:
    // ------------------------------------------------------------------------------    
    public void initializeParticles() {
        desiredKineticEnergy = initialKineticEnergy;
        setRandomPositions();
        setVelocities();
        computeAccelerations();
        paintAtomColors();
    }

    // ------------------------------------------------------------------------------
    // Place particles on a rectangular lattice:
    // ------------------------------------------------------------------------------    
    public void setRectangularLattice() {
        int nx = (int) Math.sqrt(N);
        int ny = nx;
        double dx = boxWidth / nx;
        double dy = boxWidth / ny;
        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy < ny; iy++) {
                int i = ix + iy * ny;
                x[i] = dx * (ix + 0.5);
                y[i] = dy * (iy + 0.5);
            }
        }
    }

    // ------------------------------------------------------------------------------
    // Place particles on a triangular lattice:
    // ------------------------------------------------------------------------------    
    public void setTriangularLattice() {
        int nx = 25;             //(int)Math.sqrt(N);
        int ny = 25;             //nx;
        double dx = 1.1;         //boxWidth/nx;
        double dy = 1.1 * 0.866; //boxWidth/ny;
        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy < ny; iy++) {
                int i = iy + ix * ny;
                y[i] = dy * (iy + 0.5) + 12.0;
                if (iy % 2 == 0) {
                    x[i] = dx * (ix + 0.25) + 12.0;
                } else {
                    x[i] = dx * (ix + 0.75) + 12.0;
                }
            }
        }
    }

    // ------------------------------------------------------------------------------
    // Place particles in a random distribution:
    // ------------------------------------------------------------------------------    
    public void setRandomPositions() {
        for (int i = 0; i < N; i++) {
            x[i] = boxWidth * (Math.random());
            y[i] = boxWidth * (Math.random());
            for (int j = 0; j < i; j++) {
                double dx = x[i] - x[j];
                double dy = y[i] - y[j];
                double rijSq = dx * dx + dy * dy;
                if (rijSq < 1.0) {
                    i = i - 1;
                }
            }
        }

    }

    // ------------------------------------------------------------------------------
    // Set the initial velocities of the atoms:
    // ------------------------------------------------------------------------------    
    public void setVelocities() {

        double vxSum = 0.0;
        double vySum = 0.0;

        // assign random initial velocities:
        for (int i = 0; i < N; i++) {
            vx[i] = Math.random() - 0.5;
            vy[i] = Math.random() - 0.5;
            vxSum += vx[i];
            vySum += vy[i];
        }

        // zero center of mass momentum:
        double vxcm = vxSum / N;
        double vycm = vySum / N;
        for (int i = 0; i < N; i++) {
            vx[i] -= vxcm;
            vy[i] -= vycm;
        }

        // rescale velocities to get desired initial kinetic energy:
        double v2sum = 0.0;
        for (int i = 0; i < N; i++) {
            v2sum += vx[i] * vx[i] + vy[i] * vy[i];
        }
        double kineticEnergyPerParticle = 0.5 * v2sum / N;
        double rescale = Math.sqrt(desiredKineticEnergy / kineticEnergyPerParticle);
        for (int i = 0; i < N; i++) {
            vx[i] *= rescale;
            vy[i] *= rescale;
        }
    }

    // ------------------------------------------------------------------------------
    // Check if system needs to be re-initialized:
    // ------------------------------------------------------------------------------    
    synchronized void checkControls() {
        if (reinit) {
            reinit = false;
            initializeParticles();
            theCanvas.repaint();
            gravitySwitch.setSelected(true);
            slowmotionSwitch.setSelected(false);
        }
    }

    // ------------------------------------------------------------------------------
    // Here's the method that actually runs the simulation:
    // ------------------------------------------------------------------------------     
    @Override
    public void run() {
        int stepCount = 0;

        while (true) {

            int sleepDuration = 20; //1;                
            int maxRealTime = 20; //200;
            double maxSimTime = maxRealTime / 50.0;

            if (slowmotionSwitch.isSelected()) {
                sleepDuration = 40;
                maxRealTime = 10;
                maxSimTime = maxRealTime / 50.0;
            }

            // check controls (i.e. buttons, dials, etc):
            checkControls();

            // Now for the action:
            if (running) {
                // time & step calculations
                long realTimeLimit = System.currentTimeMillis() + maxRealTime;
                double simTimeLimit = t + maxSimTime;

                // update simulation  !!!!
                while ((t < simTimeLimit) && (System.currentTimeMillis() < realTimeLimit)) {
                    singleStep();
                    stepCount++;
                }

                // draw the atoms on the screen
                paintAtomColors();
                theCanvas.repaint();

                // sleep a little
                try {
                    Thread.sleep(sleepDuration);
                } catch (InterruptedException e) {
                }

                // write to console which class is being used 
                writeMessagesToConsole();

            }
        }
    }

    // ------------------------------------------------------------------------------
    // Execute one time step using the Verlet algorithm (from Gould and Tobochnik):
    // ------------------------------------------------------------------------------    
    synchronized void singleStep() {
        // velocity-verlet algorithm:
        updatePositions();
        updateVelocities();
        computeAccelerations();
        updateVelocities();

        // update time
        t += dt;

        // rescale velocities to get desired 
        // kinetic energy:
        rescaleVelocities();

    }

    // ------------------------------------------------------------------------------
    // Update atom positions:
    // ------------------------------------------------------------------------------
    void updatePositions() {
        double dtSquaredOver2 = 0.5 * dt * dt;
        for (int i = 0; i < N; i++) {
            x[i] += vx[i] * dt + ax[i] * dtSquaredOver2;
            y[i] += vy[i] * dt + ay[i] * dtSquaredOver2;
        }
        derivedClassPos = false;
    }

    // ------------------------------------------------------------------------------
    // Execute one time step using the Verlet algorithm (from Gould and Tobochnik):
    // ------------------------------------------------------------------------------
    void updateVelocities() {
        double dtOver2 = 0.5 * dt;
        for (int i = 0; i < N; i++) {
            vx[i] += ax[i] * dtOver2;  // finish updating velocity with new acceleration
            vy[i] += ay[i] * dtOver2;
        }
        derivedClassVel = false;
    }

    // ------------------------------------------------------------------------------
    // Update atom positions:
    // ------------------------------------------------------------------------------
    void rescaleVelocities() {

        // sum the velocities squared:
        double v2sum = 0.0;
        for (int i = 0; i < N; i++) {
            v2sum += vx[i] * vx[i] + vy[i] * vy[i];
        }

        // rescale velocities based on ratio of desired/actual K.E.:
        double kineticEnergyPerParticle = 0.5 * v2sum / N;
        double rescale = Math.sqrt(desiredKineticEnergy / kineticEnergyPerParticle);
        for (int i = 0; i < N; i++) {
            vx[i] *= rescale;
            vy[i] *= rescale;
        }

        // update the thermostat bar:
        thermostat.setValue((int) (kineticEnergyPerParticle / maxKineticEnergy * 100));
    }

    // ------------------------------------------------------------------------------
    // Determine an atom color based on its kinetic energy:
    // ------------------------------------------------------------------------------ 
    void paintAtomColors() {
        // loop over atoms...
        for (int i = 0; i < N; i++) {
            mColor[i] = atomSpeedColor(i);
        }
    }

    // ------------------------------------------------------------------------------
    // Determine an atom color based on its kinetic energy
    // (this is my color spectrum):
    // ------------------------------------------------------------------------------ 
    Color atomSpeedColor(int i) {
        double iKE = 0.5 * (vx[i] * vx[i] + vy[i] * vy[i]) / maxKineticEnergy;
        iKE = Math.sqrt(iKE);
        if (iKE > 1.0) {
            iKE = 1.0;
        }
        iKE = 0.65f + 0.35f * iKE;
        return Color.getHSBColor((float) iKE, 1.0f, 1.0f);
    }

    // ------------------------------------------------------------------------------
    // Determine an atom color based on its kinetic energy
    // (this is Dan Schroeder's spectrum):
    // ------------------------------------------------------------------------------
    Color speedColor(int i) {
        int colorCount = speedColorList.length;
        double speed = Math.sqrt(vx[i] * vx[i] + vy[i] * vy[i]);
        final double speedLimit = 3.0;     // above this (squared) speed, all get the same color
        if (speed >= speedLimit) {
            return speedColorList[colorCount - 1];
        } else {
            return speedColorList[(int) (speed * colorCount / speedLimit)];
        }
    }

    // ------------------------------------------------------------------------------
    // Compute accelerations of all atoms from current positions:
    // ------------------------------------------------------------------------------ 
    void computeAccelerations() {
        // first check for bounces off walls, and include gravity (if any):
        computeWallForces();

        // Now compute interaction forces (Lennard-Jones potential).
        // This is where the program spends most of its time (when N is reasonably large),
        // so we carefully avoid unnecessary calculations and array lookups.
        computeAtomPairForces();

    }

    // ------------------------------------------------------------------------------
    // Compute wall forces exerted on atoms (wall forces are modeled as an 
    // elastic surface):
    // ------------------------------------------------------------------------------ 
    void computeWallForces() {
        for (int i = 0; i < N; i++) {
            // wall forces...
            if (x[i] < 0.5) {
                ax[i] = wallStiffness * (0.5 - x[i]);
            } else if (x[i] > (boxWidth - 0.5)) {
                ax[i] = wallStiffness * (boxWidth - 0.5 - x[i]);
            } else {
                ax[i] = 0.0;
            }
            if (y[i] < 0.5) {
                ay[i] = (wallStiffness * (0.5 - y[i]));
            } else if (y[i] > (boxWidth - 0.5)) {
                ay[i] = (wallStiffness * (boxWidth - 0.5 - y[i]));
            } else {
                ay[i] = 0;
            }

            // gravitational force...
            if (gravitySwitch.isSelected()) {
                ay[i] -= gravity;
            }
        }
    }

    // ------------------------------------------------------------------------------
    // Compute atom pair interactions based on the user-defined function:
    // ------------------------------------------------------------------------------ 
    void computeAtomPairForces() {
        // leave this for the user to define.... :)
        double dx, dy;  // separations in x and y directions
        double dx2, dy2, rSquared, rSquaredInv, attract, repel, fOverR, fx, fy;

        for (int i = 1; i < N; i++) {
            for (int j = 0; j < i; j++) {        // loop over all distinct pairs
                dx = x[i] - x[j];
                dx2 = dx * dx;
                if (dx2 < forceCutoff2) {        // make sure they're close enough to bother
                    dy = y[i] - y[j];
                    dy2 = dy * dy;
                    if (dy2 < forceCutoff2) {
                        rSquared = dx2 + dy2;
                        computeLennardJones(i, j, rSquared);
                    }
                }
            }
        }
        derivedClassSep = false;
    }

    // ------------------------------------------------------------------------------
    // Compute Lennard-Jones forces for the i-j pair:
    // ------------------------------------------------------------------------------ 
    void computeLennardJones(int i, int j, double rSquared) {
        double rSquaredInv, attract, repel, fOverR, fx, fy;
        double dx, dy;

        if (rSquared < forceCutoff2) {
            rSquaredInv = 1.0 / rSquared;
            attract = rSquaredInv * rSquaredInv * rSquaredInv;
            repel = attract * attract;
            fOverR = 24.0 * ((2.0 * repel) - attract) * rSquaredInv;
            dx = x[i] - x[j];   // note: I re-calculate dx and dy here
            dy = y[i] - y[j];   //       It is not the most efficient, but good for making a separate subroutine.
            fx = fOverR * dx;
            fy = fOverR * dy;
            ax[i] += fx;  // add this force on to i's acceleration (mass = 1)
            ay[i] += fy;
            ax[j] -= fx;  // Newton's 3rd law
            ay[j] -= fy;
        }
    }

    // ------------------------------------------------------------------------------
    // Let the user know which class is being used for each function:
    // ------------------------------------------------------------------------------ 
    void writeMessagesToConsole() {
        if (firstPrintOff == true) {
            // tell user which class 'updatePositions' is running from:
            if (derivedClassPos == true) {
                System.out.println("updatePositions: Running your class!!!!  Great Job, Java Superstar!!!!");
            } else {
                System.out.println("updatePositions: Running base class");
            }

            // tell user which class 'updateVelocities' is running from:
            if (derivedClassVel == true) {
                System.out.println("updateVelocities: Running your class!!!!  Great Job, Java Superstar!!!!");
            } else {
                System.out.println("updateVelocities: Running base class");
            }

            // tell user which class 'computeAtomPairForces' is running from:
            if (derivedClassSep == true) {
                System.out.println("computeAtomPairForces: Running your class!!!!  Great Job, Java Superstar!!!!");
            } else {
                System.out.println("computeAtomPairForces: Running base class");
            }

        }

        firstPrintOff = false;
    }

// ------------------------------------------------------------------------------
// This inner class provides the drawing space and and does all the drawing:
// ------------------------------------------------------------------------------ 
    class MDCanvas extends Canvas {

        // off-screen image for double-buffering, to prevent flicker:
        Image offScreenImage;
        Graphics offScreenGraphics;

        Cursor handCursor = new Cursor(Cursor.HAND_CURSOR);
        Cursor crosshairCursor = new Cursor(Cursor.CROSSHAIR_CURSOR);
        Font bigFont = new Font(null, Font.BOLD, 20);
        int circleSizeCorrection = 0;   // see constructor and drawOffScreenImage method
        boolean dragging = false;       // true when a molecule is being dragged
        boolean drawingBond = false;    // true when we're drawing a new bond
        int newBondIndex;               // index of molecule from which we're drawing a bond
        int bondEndX, bondEndY;         // pixel coordinates of mouse location while drawing a new bond

        // carefully constructed list of colors for indicating speeds:
        final Color[] speedColorList = {new Color(20, 0, 180), new Color(60, 0, 170), new Color(80, 0, 160), new Color(100, 0, 150),
            new Color(120, 0, 120), new Color(140, 0, 80), new Color(160, 0, 0), new Color(180, 0, 0),
            new Color(200, 0, 0), new Color(230, 0, 0), Color.RED, new Color(255, 60, 0),
            new Color(255, 90, 0), new Color(255, 120, 0), new Color(255, 150, 0), new Color(255, 180, 0),
            new Color(255, 210, 0), new Color(255, 230, 0), Color.YELLOW, new Color(255, 255, 120)};

        // canvas constructor method:
        MDCanvas() {
            setSize(canvasWidth, canvasWidth);
            setCursor(handCursor);
            if (System.getProperty("os.name").equals("Mac OS X") && (System.getProperty("os.version").charAt(3) < '5') && (System.getProperty("java.vm.version").charAt(2) >= '4')) {
                circleSizeCorrection = -1; // Java 1.4+ on Mac OS X 10.4- draws circles one pixel too big!
            }
        }

        // draw the off-screen image, for later copying to screen:
        // (This is where we check the popup menus for the selected colors.  Note that the last
        // two entries in the molecule color menu are hard-coded for "random" and "by speed".)
        void drawOffScreenImage(Graphics g) {
            g.setColor(Color.BLACK);
            g.fillRect(0, 0, canvasWidth, canvasWidth);
            int circleSize = pixelsPerUnit + circleSizeCorrection;
            for (int i = 0; i < N; i++) {
                int screenx = xToScreen(x[i] - 0.5);     // circle is drawn from upper-left corner
                int screeny = yToScreen(y[i] + 0.5);
                g.setColor(mColor[i]);
                if (pixelsPerUnit < 5) {
                    g.fillRect(screenx, screeny, pixelsPerUnit, pixelsPerUnit);
                } else {
                    g.fillOval(screenx, screeny, circleSize, circleSize);
                }
            }
        }

        // return a color to indicate the speed of molecule i:
        Color speedColor(int i) {
            int colorCount = speedColorList.length;
            double speed = Math.sqrt(vx[i] * vx[i] + vy[i] * vy[i]);
            final double speedLimit = 3.0;     // above this (squared) speed, all get the same color
            if (speed >= speedLimit) {
                return speedColorList[colorCount - 1];
            } else {
                return speedColorList[(int) (speed * colorCount / speedLimit)];
            }
        }

        // override update to skip painting background (improves performance and reduces flicker)
        @Override
        public void update(Graphics g) {
            paint(g);
        }

        // paint method draws the off-screen image first, then blasts it to screen:
        @Override
        public void paint(Graphics g) {
            if (offScreenImage == null) {
                offScreenImage = createImage(canvasWidth, canvasWidth);    // can't do this until paint is first called
                offScreenGraphics = offScreenImage.getGraphics();
            }
            drawOffScreenImage(offScreenGraphics);
            g.drawImage(offScreenImage, 0, 0, null);
        }

        // Look for a molecule at given pixel coordinates and return its index, or -1 if not found:
        int findMolecule(int pixelX, int pixelY) {
            double mX = xToActual(pixelX);      // convert coordinates to simulation units
            double mY = yToActual(pixelY);
            final double radius2 = 0.55 * 0.55; // square of molecular radius (with a little extra)
            boolean found = false;              // no molecule found yet
            int i = 0;                          // molecule index
            while (i < N) {
                double dx = mX - x[i];
                double dx2 = dx * dx;
                if (dx2 < radius2) {
                    double dy = mY - y[i];
                    double dy2 = dy * dy;
                    if (dy2 < radius2) {
                        double r2 = dx2 + dy2;
                        if (r2 < radius2) {
                            found = true;
                            break;
                        }
                    }
                }
                i++;
            }
            if (found) {
                return i;
            } else {
                return -1;
            }
        }

        // conversions from simulation coordinates to screen coordinates:
        int xToScreen(double xActual) {
            return (int) Math.round(xActual * pixelsPerUnit);
        }

        int yToScreen(double yActual) {
            return canvasWidth - (int) Math.round(yActual * pixelsPerUnit);
        }

        // conversions from screen coordinates to simulation coordinates:
        double xToActual(int xScreen) {
            return xScreen * 1.0 / pixelsPerUnit;
        }

        double yToActual(int yScreen) {
            return (canvasWidth - yScreen) * 1.0 / pixelsPerUnit;
        }

        // not needed since we handle press & release separately
        public void mouseClicked(MouseEvent e) { }  

    } // end of inner class MDCanvas

}
