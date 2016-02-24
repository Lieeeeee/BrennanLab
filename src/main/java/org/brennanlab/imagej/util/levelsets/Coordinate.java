package org.brennanlab.imagej.util.levelsets;

/**
 * From fiji
 *
 * @author fiji
 */
public class Coordinate {
    private int x; // column
    private int y; // row

    public Coordinate(int x, int y) {
        this.x = x;
        this.y = y;
    }

    public int getX() {
        return x;
    }

    public void setX(int x) {
        this.x = x;
    }

    public int getY() {
        return y;
    }

    public void setY(int y) {
        this.y = y;
    }

    public boolean equals(Coordinate p) {
        return (x == p.getX() && y == p.getY());
    }

    public double distance(Coordinate p) {
        return Math.sqrt((x - p.getX()) * (x - p.getX()) + (y - p.getY())
                * (y - p.getY()));
    }
}