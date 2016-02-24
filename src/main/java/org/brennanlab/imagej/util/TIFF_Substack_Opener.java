package org.brennanlab.imagej.util;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Menus;
import ij.gui.GenericDialog;
import ij.io.*;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;

public class TIFF_Substack_Opener implements PlugIn {
    private String userInput, stackTitle, num, strA, strB;
    private int rangeStart, rangeEnd, range, currSlice, count, width, height;
    private String path, name, directory;
    private int stackSize;
    private boolean bAbort;
    private boolean skipinput = false;
    private int[] numList;
    private FileInfo[] fi;
    private Integer obj;
    private ImagePlus impSubstack;

    public TIFF_Substack_Opener() {

    }

    public TIFF_Substack_Opener(String directory, String name) {
        this.name = name;
        this.directory = directory;
        this.path = directory + name;
        this.skipinput = true;
        TiffDecoder td = new TiffDecoder(directory, name);
        try {
            FileInfo[] fi = td.getTiffInfo();
            width = fi[0].width;
            height = fi[0].height;
            stackSize = fi[0].nImages;
        } catch (IOException e) {
            bAbort = true;
            return;
        }
    }

    InputStream createInputStream(FileInfo fi) throws IOException,
            MalformedURLException {
        if (fi.inputStream != null)
            return fi.inputStream;
        else if (fi.url != null && !fi.url.equals(""))
            return new URL(fi.url + fi.fileName).openStream();
        else {
            File f = new File(fi.directory + fi.fileName);
            if (f == null || f.isDirectory())
                return null;
            else {
                InputStream is = new FileInputStream(f);
                if (fi.compression >= FileInfo.LZW)
                    is = new RandomAccessStream(is);
                return is;
            }
        }
    }

    void getInput() {

        GenericDialog gd = new GenericDialog("Substack Opener", IJ
                .getInstance());
        gd.addMessage("There are a total of " + stackSize
                + " totalSlices in the stack you have selected");
        gd
                .addMessage("Enter either range (e.g. 2-14) or a list (e.g. 7,9,25,27)");
        gd.addStringField("totalSlices", "1-200", 50);
        gd.showDialog();
        if (gd.wasCanceled()) {
            bAbort = true;
            return;
        }
        userInput = gd.getNextString();
        if (userInput.length() == 0) {
            IJ.error("Input required.");
            bAbort = true;
            return;
        }
    }

    public ImagePlus returnRange(int start, int stop) {
        this.rangeStart = start;
        this.rangeEnd = stop;
        range = rangeEnd - rangeStart + 1;
        return stackRange(rangeStart, "");
    }

    public ImagePlus returnResult() {
        return impSubstack;
    }

    public void run(String arg) {

        if (arg.length() > 0) {
            String[] args = arg.split(";");
            for (int i = 0; i < args.length; i++) {
                String[] param = args[i].split("=");
                if (param[0].equals("directory"))
                    directory = param[1];
                else if (param[0].equals("userInput"))
                    userInput = param[1];
                else if (param[0].equals("name"))
                    name = param[1];

            }
            if (name.length() > 0 && directory.length() > 0
                    && userInput.length() > 0) {
                skipinput = true;
            }
        }

        // First ask the user to choose a file, then ask the user to choose
        // totalSlices to open
        if (skipinput) {
            if (name != null) {
                path = directory + name;
                boolean error = false;
                if (!error)
                    Menus.addOpenRecentItem(path);
            }

            TiffDecoder td = new TiffDecoder(directory, name);
            try {
                fi = td.getTiffInfo();
                width = fi[0].width;
                height = fi[0].height;
                stackSize = fi[0].nImages;
            } catch (IOException e) {
                bAbort = true;
                return;
            }

        } else {
            selectFile();
            if (name == null || directory == null)
                return;
            if (bAbort)
                return;
            // set the width, height, number of totalSlices

            getInput();
            if (userInput == null)
                return;
            bAbort = false;

            if (bAbort)
                return;
        }
        stackTitle = "Substack (" + userInput + ")";
        if (stackTitle.length() > 25) {
            int idxA = stackTitle.indexOf(",", 18);
            int idxB = stackTitle.lastIndexOf(",");
            if (idxA >= 1 && idxB >= 1) {
                strA = stackTitle.substring(0, idxA);
                strB = stackTitle.substring(idxB + 1);
                stackTitle = strA + ", ... " + strB;
            }
        }

        try {
            int idx1 = userInput.indexOf("-");
            if (idx1 >= 1) { // input displayed in range
                String rngStart = userInput.substring(0, idx1);
                String rngEnd = userInput.substring(idx1 + 1);
                obj = new Integer(rngStart);
                rangeStart = obj.intValue();
                obj = new Integer(rngEnd);
                rangeEnd = obj.intValue();
                range = rangeEnd - rangeStart + 1;
                impSubstack = stackRange(rangeStart, stackTitle);
            } else {
                count = 1; // count # of totalSlices to extract
                for (int j = 0; j < userInput.length(); j++) {
                    char ch = Character.toLowerCase(userInput.charAt(j));
                    if (ch == ',') {
                        count += 1;
                    }
                }

                numList = new int[count];
                for (int i = 0; i < count; i++) {
                    int idx2 = userInput.indexOf(",");
                    if (idx2 > 0) {
                        num = userInput.substring(0, idx2);
                        obj = new Integer(num);
                        numList[i] = obj.intValue();
                        userInput = userInput.substring(idx2 + 1);
                    } else {
                        num = userInput;
                        obj = new Integer(num);
                        numList[i] = obj.intValue();
                    }
                }
                impSubstack = stackList(numList, stackTitle);
            }
        } catch (NumberFormatException e) {
            IJ.error("Improper input:\n" + userInput);
            return;
        }
        this.showResult();
    }

    void selectFile() {
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
            width = fi[0].width;
            height = fi[0].height;
            stackSize = fi[0].nImages;
        } catch (IOException e) {
            bAbort = true;
            return;
        }

    }

    public void showResult() {
        impSubstack.show();
    }

    ImagePlus stackList(int[] numList, String stackTitle) { // extract specific
        // totalSlices
        ImagePlus impSub = null;
        try {
            Opener op = new Opener();
            currSlice = numList[0];
            ImageProcessor ip2 = op.openTiff(path, currSlice).getProcessor();
            ImageStack stack2 = new ImageStack(width, height, ip2
                    .getColorModel());
            stack2.addSlice(stackTitle + " " + currSlice, ip2);

            for (int i = 1; i < count; i++) {
                currSlice = numList[i];
                ImageProcessor ip3 = op.openTiff(path, currSlice)
                        .getProcessor();
                ip3 = ip3.crop();
                stack2.addSlice(stackTitle + " " + currSlice, ip3);
                IJ.showProgress(i, numList.length);
            }

            impSub = new ImagePlus("", stack2);

        } catch (IllegalArgumentException e) {
            IJ.error("Argument out of range: " + userInput);
        }
        return impSub;
    }

    ImagePlus stackRange(int currSlice, String stackTitle) { // extract range of
        TiffDecoder td = new TiffDecoder(directory, name);
        if (IJ.debugMode)
            td.enableDebugging();
        FileInfo[] info = null;
        try {
            info = td.getTiffInfo();
        } catch (IOException e) {
            String msg = e.getMessage();
            if (msg == null || msg.equals(""))
                msg = "" + e;
            IJ.error("TiffDecoder", msg);
            return null;
        }
        if (info == null)
            return null;
        FileInfo fi = info[0];
        if (info.length == 1 && fi.nImages > 1) {
            if (rangeStart < 1 || rangeStart > fi.nImages)
                throw new IllegalArgumentException("N out of 1-" + fi.nImages
                        + " range");
            long size = fi.width * fi.height * fi.getBytesPerPixel();
            fi.longOffset = fi.getOffset() + (rangeStart - 1)
                    * (size + fi.gapBetweenImages);
            fi.offset = 0;
            fi.nImages = range;
        } else {
            if (rangeStart < 1 || rangeEnd > info.length)
                throw new IllegalArgumentException("N out of 1-" + info.length
                        + " range");
            fi.longOffset = info[rangeStart - 1].getOffset();
            fi.offset = 0;
            fi.stripOffsets = info[rangeStart - 1].stripOffsets;
            fi.stripLengths = info[rangeStart - 1].stripLengths;
        }
        FileOpener fo = new FileOpener(fi);

        ImagePlus impSub = fo.open(false);
        for (int i = rangeStart; i <= rangeEnd; i++)
            impSub.getStack().setSliceLabel(fi.fileName + " slice " + i,
                    i - rangeStart + 1);

        impSub.setTitle(fi.fileName + " " + rangeStart + "-" + rangeEnd);
        return impSub;

    }

}
