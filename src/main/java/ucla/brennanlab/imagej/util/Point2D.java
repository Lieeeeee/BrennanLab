package ucla.brennanlab.imagej.util;

public class Point2D{
    public int x;
    public int y;
    public Point2D(int x, int y) {
        this.x = x;
        this.y = y;
    }
    public Point2D(int[] x){
        this.x = x[0];
        this.y = x[1];
    }
    public int getX() {
        return x;
    }

    public int getY() {
        return y;
    }

    public boolean equals(Object other){
        if (!(other instanceof Point2D))
            return false;
        Point2D otherPoint = (Point2D) other;
        return x==otherPoint.x && y==otherPoint.y;
    }

    public int hashCode(){
        return x ^ y;
    }

    public double distanceFrom(Point2D other){
        int dx = x-other.x;
        int dy = y-other.y;
        return Math.sqrt(dx*dx+dy*dy);
    }

    public String toString(){
        return "("+x+","+y+")";
    }

}
