package ucla.brennanlab.imagej.util;

public class Point2D {
    public int x,y;
    public Point2D(int x , int y){
        this.x = x;
        this.y = y;
    }
    public Point2D(int[] x){
        this.x = x[0];
        this.y = x[1];
    }
    public double distance(Point2D x){
        return Math.sqrt( (this.x-x.x)^2 + (this.y-x.y)^2 );
    }
}
