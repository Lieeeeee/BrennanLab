package ucla.brennanlab.imagej.util;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

class Pixel_scatterplot implements PlugIn {
	private static String title1 = "";
	private static String title2 = "";

	@Override
	public void run(String arg) {
		int[] wList = WindowManager.getIDList();
		if (wList == null) {
			IJ.noImage();
			return;
		}
		String[] titles = new String[wList.length];
		for (int i = 0; i < wList.length; i++) {
			ImagePlus imp = WindowManager.getImage(wList[i]);
			if (imp != null)
				titles[i] = imp.getTitle();
			else
				titles[i] = "";
		}
		GenericDialog gd = new GenericDialog("Image Calculator", IJ
				.getInstance());
		String defaultItem;
		if (title1.equals(""))
			defaultItem = titles[0];
		else
			defaultItem = title1;
		gd.addChoice("Image1:", titles, defaultItem);
		if (title2.equals(""))
			defaultItem = titles[0];
		else
			defaultItem = title2;
		gd.addChoice("Image2:", titles, defaultItem);
		// gd.addStringField("Result:", "Result", 10);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		int index1 = gd.getNextChoiceIndex();
		title1 = titles[index1];
		int index2 = gd.getNextChoiceIndex();
		// String resultTitle = gd.getNextString();

		title2 = titles[index2];
		@SuppressWarnings("unused")
		ImagePlus img1 = WindowManager.getImage(wList[index1]);
		@SuppressWarnings("unused")
		ImagePlus img2 = WindowManager.getImage(wList[index2]);

	}

}