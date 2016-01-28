package ucla.brennanlab;

import com.google.common.collect.Sets;
import org.delaunay.algorithm.Triangulation.InvalidVertexException;
import org.delaunay.dtfe.DtfeTriangulationMap;
import org.delaunay.dtfe.interpolation.InterpolationStrategies;
import org.delaunay.model.Vector;
import org.delaunay.model.Vertex;

import java.util.HashMap;
import java.util.LinkedHashSet;


public class DelaunayInterpolator extends DtfeTriangulationMap{
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
        InterpolationStrategies.createNaturalNeighbor();
        for(int x=0; x<width; x++){
            for(int y=0; y<height; y++){
                speeds[x][y] = (float)this.getInterpolatedDensity(new Vector(x,y),
                        InterpolationStrategies.createNaturalNeighbor());
            }
        }

        return speeds;
    }

    public double getDensity(Vertex v){
        return this.values.get(v);
    }
}