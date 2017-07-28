package game;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseEvent;
import java.awt.Image;
import java.awt.Graphics;
import java.io.File;
import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.swing.Timer;
import java.awt.Color;
import javax.swing.JPanel;
import java.util.Random;

/* Information about the sprite */
class Model {

    // sprite information
    protected int s_x;     // sprite x
    protected int s_y;     // spirte y
    protected double vel;  // sprite velocity
    protected int width;   // width of image
    protected int height;  // height of image
    protected String name; // name of file
    protected Image sprite;

    /**
     * @param time sleep the thread for time milliseconds
     */
    protected void sleep(int time) {
        try {
            Thread.sleep(time);
        } catch (Exception e) {
            System.out.println("Could not sleep.");
            e.printStackTrace(System.err);
            System.exit(1);
        }
    }

    /**
     * @return x position of the right of the sprite
    */
    protected int getRight() {
        return s_x + width;
    }

    /**
     * @return y position of the bottom of the sprite
     */
    protected int getBottom() {
        return s_y + height;
    }
    
    /**
     * Loads a sprite
     * @param name name of sprite to load
     */
    protected void loadImage(String name) {
        try {
            sprite = ImageIO.read(new File("images/" + name));
        } catch (Exception e) {
            System.out.println("Current dir: " + System.getProperty("user.dir"));
            System.out.println("Could not load " + name);
            e.printStackTrace(System.err);
            System.exit(1);
        }
    }
}

/* Bird class */
class Bird extends Model {

    private View v;            // sprite looks
    private final int gravity; // gravity

    /**
     * Initialize Bird variables
     */
    protected Bird() {
        gravity = 2;
        s_x = 50;
        s_y = Game.WINDOW_SIZE / 2;
        vel = 0;
        name = "bird1.png";
    }

    /**
     * Gets the View object so it can be used by bird
     */
    protected void getView(View v) {
        this.v = v;
        name = v.n;
        width = sprite.getWidth(v);
        height = sprite.getHeight(v);
    }

    /**
     * Paints the sprite on screen 
     */
    private void paint() {
        loadImage(this.name);
        v.repaint();
    }

    /** 
     * Updates the location of the sprite and which sprite is in use 
     */
    protected void update() {
        // if the sprite is not currently rising, fall with gravity
        s_y += (int) (vel + 0.5 * gravity);
        vel += gravity;
        paint();
    }

    /**
     * On mouse click, raise the sprite and change the image 
     */
    protected void flap() {
        // loads other sprite on click
        if (name.equals("bird2.png")) {
            name = "bird1.png";
        } else {
            name = "bird2.png";
        }

        // raises the sprite up
        vel = -14;
    }
}

/* Pipe class */
class Pipe extends Model {

    // variables
    private final int min;
    private final int max;
    private final boolean invert;

    /**
     * Initialize variables for pipe objects
     * @param min Used to randomly pick a y-location for the pipe to spawn
     * @param max Used to randomly pick a y-location for the pipe to spawn
     * @param invert The top pipe will have its random location inverted
     * @param r Random number generation 
     * @param name name of image
     */
    protected Pipe(int min, int max, boolean invert, Random r, String name) {
        // get pipe information
        this.min = min;
        this.max = max;
        this.invert = invert;
        this.name = name;
        loadImage(name);
        reset(r);

        // default information
        width = 55;
        height = 400;
        vel = 5;
    }

    /**
     * Generate an x and y position for newly generated pipes 
     * @param r Random number generation
     */
    protected void reset(Random r) {
        // get x and y positions
        s_y = r.nextInt(max - min) + min;
        s_x = Game.WINDOW_SIZE - 50;
        
        // if top pipe, then switch generated y to a negative number
        if (invert) {
            s_y *= -1;
        }
        
        // make the pipes faster
        vel += 0.2;
    }
}

/* How the sprite looks */
class View extends JPanel {

    // sprite information
    private final Model m;  // sprite information such as x and y position
    protected String n;     // name of sprite
    
    // pipe information
    private final Random r;
    private final Pipe[] p;

    /**
     * View constructor: loads the sprite and the pipes
     * @param m stores the model information so we can later access information like x and y position
     * @param name name of initial sprite to load
     */
    protected View(Model m, String name) {
        p = new Pipe[2];
        r = new Random();
        
        // sprite information
        this.m = m;
        n = name;
        m.loadImage(n);

        // pipe information: do while the gap is too small or the pipes are overlapping
        do {
            p[0] = new Pipe(0, 350, true, r, "tube_down.png"); // top pipe
            p[1] = new Pipe(350, 650, false, r, "tube_up.png");  // bottom pipe
        } while (pipesOverlapping() || pipesFarApart());
    }
    
    /**
     * TODO: Eventually move into Pipe
     * @return whether or not the pipes are overlapping
     */
    private boolean pipesOverlapping() {
        return (p[1].s_y - p[0].getBottom() < 150 || p[0].getBottom() > p[1].s_y);
    }
    
    /**
     * TODO: Eventually move into Pipe
     * @return whether or not the pipes are too far apart for a challenge 
     */
    private boolean pipesFarApart() {
        return (p[1].s_y - p[0].getBottom() > 300);
    }

    /**
     * TODO: Eventually put parts into Pipe?
     * Draws the pipes on screens, moves the pipes, and will call reset if needed
     * @param g needed to be able to draw the pipes
     */
    private void drawPipe(Graphics g) {
        // for each pipe
        for (Pipe pipe : p) {
            // draw the pipe
            g.drawImage(pipe.sprite, pipe.s_x, pipe.s_y, null);

            // move the pipe
            pipe.s_x -= pipe.vel;
        }

        // if pipe has reached edge of screen, generate some new pipes
        // pipes should be at same x position so only need to check one pipe
        if (p[0].s_x <= -100) {
            // do while the gap is too small or the pipes are overlapping
            do {
                p[0].reset(r);
                p[1].reset(r);
            } while (pipesOverlapping() || pipesFarApart());
        }
    }

    /**
     * @return whether or not user is colliding with edge of screen or pipes 
     */
    private boolean collision() {
        // leeway allows the bird to collide with the pipe a bit because of the invisible border around the bird
        final int leeway = 15;
        return ((m.getRight() >= p[0].s_x && m.s_x <= p[0].getRight()) && (m.s_y <= p[0].getBottom() - leeway || m.getBottom() - leeway >= p[1].s_y)) || (m.getBottom() >= Game.WINDOW_SIZE + 15 || m.s_y <= 0 - leeway);
    }

    /**
     * Draws to screen. Also, calls methods to check collisions
     * @param g needed to draw to screen
     */
    @Override
    public void paintComponent(Graphics g) {
        // draw the screen
        //super.paintComponent(g); // uncomment if things start to act wierd
        setBackground(Color.cyan);
        drawPipe(g);
        
        // prints lose to the screen if user has collided with pipes or edges of screens
        if (collision()) {
            System.out.println("Lose");
            m.sleep(5000);
            System.exit(0);
        }

        // draw sprite
        if (m.sprite != null && p[0].sprite != null && p[1].sprite != null) {
            g.drawImage(m.sprite, m.s_x, m.s_y, null);
            drawPipe(g);
        } else {
            throw new RuntimeException("One of the sprites is null\n" + m.sprite + "\n" + p[0].sprite + "\n" + p[1].sprite);
        }

    }
}

/* gets actions for the sprite */
class Controller implements ActionListener, MouseListener {

    private final Bird b;

    /**
     * Controller constructor
     * @param b takes in a bird object so controller can tell bird to flap on mouse click
     */
    protected Controller(Bird b) {
        this.b = b;
    }

    /**
     * If mouse is pressed, sprite flaps wings
     * @param e detects mouse presses
     */
    @Override
    public void mousePressed(MouseEvent e) {
        b.flap();
    }

    @Override
    public void mouseReleased(MouseEvent e) {
    }

    @Override
    public void mouseEntered(MouseEvent e) {
    }

    @Override
    public void mouseExited(MouseEvent e) {
    }

    @Override
    public void mouseClicked(MouseEvent e) {
    }
    
    @Override
    public void actionPerformed(ActionEvent e) {
    }
}

/* main class */
class Game extends JFrame implements ActionListener {

    // constants
    public static final int WINDOW_SIZE = 700;
    private final int TIME_INTERVAL = 40;

    // sprite information
    static View vb;
    Bird b;
    Model m;

    /** 
     * Initialize information about the game 
     */
    private Game() {
        // bird stuff
        b = new Bird();
        Controller c = new Controller(b);
        vb = new View(b, "bird1.png");
        b.getView(vb);
        vb.addMouseListener(c);
        new Timer(TIME_INTERVAL, this).start();
    }

    /**
     * Update the bird's location and repaint the screen
     * @param evt an event that takes place
     */
    @Override
    public void actionPerformed(ActionEvent evt) {
        b.update();
        repaint();
    }

    /**
     * Main method
     * @param args no command line arguments needed
     */
    public static void main(String[] args) {
        Game g = new Game();
        g.setTitle("Not Flappy Bird");
        g.setSize(WINDOW_SIZE, WINDOW_SIZE);
        g.setBackground(Color.cyan);
        g.add(vb);
        g.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        g.setVisible(true);
    }
}
