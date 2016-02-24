package org.brennanlab.imagej.csd;

import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.Wand;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

import java.awt.*;

public class SignedDistancetoRoi implements PlugInFilter {

    ImagePlus imp;

    @Override
    public void run(ImageProcessor ip) {

        Wand w = new Wand(ip);
        int minpt[] = {0, 0};
        for (int x = 0; x < ip.getWidth(); x++) {
            for (int y = 0; y < ip.getHeight(); y++) {
                if (ip.getPixelValue(x, y) < ip.getPixel(minpt[0], minpt[1])) {
                    minpt[0] = x;
                    minpt[1] = y;
                }
            }
        }
        w.autoOutline(minpt[0], minpt[1], ip.getPixel(minpt[0], minpt[1]), 0);
        Polygon p = new Polygon(w.xpoints, w.ypoints, w.npoints);
        PolygonRoi r = new PolygonRoi(p, Roi.POLYLINE);
        r.fitSpline();
        ip.setRoi(r);
        imp.setRoi(r);
        Overlay o = new Overlay(r);
        imp.setOverlay(o);
        imp.updateAndDraw();
    }

    @Override
    public int setup(String arg, ImagePlus imp) {
        // TODO Auto-generated method stub
        this.imp = imp;
        return DOES_32;
    }
}