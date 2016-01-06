package ucla.brennanlab.imagej.util;

import ij.IJ;
import ij.Menus;
import ij.gui.GenericDialog;
import ij.io.FileInfo;
import ij.io.OpenDialog;
import ij.io.TiffDecoder;
import ij.plugin.PlugIn;

import java.io.IOException;

public class Multiple_Substack_Opener implements PlugIn {

    private String userInput;
    private String path, name, directory;
    private int stackSize;

    public Multiple_Substack_Opener() {
        return;
    }

    public void run(String arg) {

        // ***************** Get file ***************************

        OpenDialog od = new OpenDialog("Open a TIFF file", "");
        directory = od.getDirectory();
        name = od.getFileName();
        if (name != null) {
            path = directory + name;
            boolean error = false;
            if (!error)
                Menus.addOpenRecentItem(path);
        }

        TiffDecoder td = new TiffDecoder(directory, name);
        try {
            FileInfo[] fi = td.getTiffInfo();

            stackSize = fi[0].nImages;
        } catch (IOException e) {
            return;
        }

        // ************** Get slice ranges ********************
        GenericDialog gd = new GenericDialog("Substack Opener", IJ
                .getInstance());
        gd.addMessage("There are a total of " + stackSize
                + " totalSlices in the stack you have selected");
        gd
                .addMessage("This plugin is very non-robust to bad user input, so please format your input correctly");
        gd
                .addMessage("Enter slice ranges separated by commas (e.g. 1-100,201-300,501-700)");
        gd.addStringField("totalSlices", 1 + "-" + stackSize, 50);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return;
        }
        userInput = gd.getNextString();
        if (userInput.length() == 0) {
            IJ.error("Input required.");
            return;
        }

        // ************* Parse input *****************

        // ************ Call TIFF_Substack_Opener******
        String[] sliceRanges = userInput.split(",");
        for (int i = 0; i < sliceRanges.length; i++) {
            String args = "name=" + name + ";directory=" + directory
                    + ";userInput=" + sliceRanges[i];
            IJ.runPlugIn("ucla.brennanlab.imagej.util.TIFF_Substack_Opener",
                    args);
        }
        return;

    }

}