package grid;

import java.util.Random;

// ball class
class Ball {
    int x;
    int y;
    int v_x;
    int v_y;
    
    /* ball constructor */
    Ball(int width, int height) {
        Random r = new Random();
        
        x = r.nextInt(width);
        y = r.nextInt(height);
        v_x = (r.nextInt() % 2 == 0) ? -1 : 1;
        v_y = (r.nextInt() % 2 == 0) ? -1 : 1;
    }
    
    /* ball constructor */
    Ball(int x, int y, int v_x, int v_y) {
        this.x = x;
        this.y = y;
        this.v_x = v_x;
        this.v_y = v_y;
    }
    
    /* move the ball */
    void move() {
        x += v_x;
        y += v_y;
    }
    
    /* bounce ball if touching edge of grid object */
    void bounce(Grid g) {
        if (x <= 0 || x >= g.width - 1) {
            v_x *= -1;
        }
        if (y <= 0 || y >= g.height - 1) {
            v_y *= -1;
        }
    }
}

// grid class
class Grid {
    int width;
    int height;
    
    /* main function */
    public static void main(String[] args) {
        // constants
        final int SEPERATOR = 15;
        final int SIZE = 10;
        
        // create the objects
        Grid grid = new Grid();
        Ball[] balls = new Ball[SIZE];
        for (int i = 0; i < SIZE; i++) {
            balls[i] = new Ball(grid.width, grid.height);
        }
        
        // game loop
        while(true) {
            // print SEPERATOR blank lines between grids
            for (int i = 0; i < SEPERATOR; i++) {
                System.out.println();
            }
            
            // draw ball
            grid.draw(balls, SIZE);
            sleep(500);
            for (int i = 0; i < SIZE; i++) {
                balls[i].bounce(grid);
                balls[i].move();
            }
        }
    }
    
    /* grid constructor */
    Grid() {
        width = 60;
        height = 20;
    }
    
    /* grid constructor */
    Grid(int width, int height) {
        this.width = width;
        this.height = height;
    }
    
    /* draws the grid with the ball objects
     * balls do not interact with each other */
    void draw(Ball[] balls, final int SIZE) {
        // prints out a width x height grid of #
        // except for a moving hole
        int count = 0;
        for (int i = 0; i < height; i++) {
            // print the grid
            for (int j = 0; j < width; j++) {
                for (int k = 0; k < SIZE; k++) {
                    // draw a ball if needed
                    if (j == balls[k].x && i == balls[k].y) {
                        System.out.print(" ");
                        count++;
                    }
                }
                
                // draw a # if needed
                if (count < width) {
                    System.out.print("#");
                    count++;
                }
            }
            
            // start the next row
            System.out.println();
            count = 0;
        }
    }
    
    /* sleep for time seconds */
    public static void sleep(int time) {
        try {
            Thread.sleep(time);
        } catch(Exception e) {
            e.printStackTrace(System.err);
            System.exit(1);
        }
    }
    
    
}
