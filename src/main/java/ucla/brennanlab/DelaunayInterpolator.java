package ucla.brennanlab;

import com.google.common.collect.Sets;
import org.delaunay.algorithm.Triangulation;
import org.delaunay.model.Triangle;
import org.delaunay.model.Vector;
import org.delaunay.model.Vertex;

import java.util.HashMap;
import java.util.LinkedHashSet;


public class DelaunayInterpolator extends Triangulation{
    private LinkedHashSet<Vertex> inputVertices = Sets.newLinkedHashSet();

    private HashMap<Vertex,Double> values = new HashMap<Vertex,Double>();

    public Vertex addVertex(double x, double y, double z) {
        Vertex vertex = new Vertex(x, y);
        this.inputVertices.add(vertex);
        this.values.put(vertex,z);
        return vertex;
    }

    public void triangulate() throws InvalidVertexException {
        super.triangulate();
    }

    public float[][] getInterpolation(int width, int height){
        float[][] speeds = new float[width][height];
        for(int x=0; x<width; x++){
            for(int y=0; y<height; y++){
                Triangle T = this.locate(new Vector(x,y));
                if(T==null){
                    speeds[x][y] = 0;
                }else {
                    // Inverse distance weighting??
                    float aval = (this.values.get(T.a)).floatValue();
                    float bval = (this.values.get(T.b)).floatValue();
                    float cval = (this.values.get(T.c)).floatValue();
                    speeds[x][y] = (this.values.get(T.a)).floatValue();
                }
            }
        }

        return speeds;
    }

    public double getDensity(Vertex v){
        return this.values.get(v);
    }
}