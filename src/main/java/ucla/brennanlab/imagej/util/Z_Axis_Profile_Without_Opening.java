package ucla.brennanlab.imagej.util;

import ij.IJ;
import ij.Menus;
import ij.io.FileInfo;
import ij.io.ImageReader;
import ij.io.OpenDialog;
import ij.io.RandomAccessStream;
import ij.io.TiffDecoder;
import ij.plugin.PlugIn;
import info.monitorenter.gui.chart.Chart2D;
import info.monitorenter.gui.chart.IAxis;
import info.monitorenter.gui.chart.IPointPainter;
import info.monitorenter.gui.chart.ITrace2D;
import info.monitorenter.gui.chart.controls.LayoutFactory;
import info.monitorenter.gui.chart.events.Chart2DActionSaveImageSingleton;
import info.monitorenter.gui.chart.pointpainters.PointPainterDisc;
import info.monitorenter.gui.chart.rangepolicies.RangePolicyFixedViewport;
import info.monitorenter.gui.chart.traces.Trace2DLtd;
import info.monitorenter.gui.chart.traces.Trace2DSimple;
import info.monitorenter.gui.chart.views.ChartPanel;
import info.monitorenter.util.Range;

import java.awt.Color;
import java.awt.Container;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Set;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.UIManager;

import ucla.brennanlab.util.jchart2d.ZoomableChartWithPointSelection;

public class Z_Axis_Profile_Without_Opening implements PlugIn, ActionListener {
	class ZoomAllAdapter implements ActionListener {
		/** The zoomable chart to act upon. */
		private ZoomableChartWithPointSelection m_zoomableChart;

		/**
		 * Creates an instance that will reset zooming on the given zoomable
		 * chart upon the triggered action.
		 * <p>
		 * 
		 * @param chart
		 *            the target to reset zooming on.
		 */
		public ZoomAllAdapter(final ZoomableChartWithPointSelection chart) {
			this.m_zoomableChart = chart;
		}

		/**
		 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
		 */
		public void actionPerformed(final ActionEvent event) {
			this.m_zoomableChart.zoomAll();
		}
	}

	public static int sum(int[] t) {
		int sum = 0;
		for (int j = 0; j < t.length; j++) {
			sum += t[j];
		}
		return sum;
	}

	public String path, name, directory;
	public int stackSize = 0;
	FileInfo[] fi;
	public JFrame frame;
	JButton snapshotImage;
	public ZoomableChartWithPointSelection chart;
	public double[] means;
	public ITrace2D trace;
	public ITrace2D filtered;

	public boolean selector;

	public Z_Axis_Profile_Without_Opening() {
		this.selector = false;
	}

	public Z_Axis_Profile_Without_Opening(boolean selector) {
		this.selector = selector;
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		// TODO Auto-generated method stub
		Object b = e.getSource();
		if (b == snapshotImage) {
			Chart2DActionSaveImageSingleton.getInstance(chart, "Save image");
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

	public void hideChart() {
		frame.setVisible(false);
	}

	public void plotData() {
		long skip = fi[0].getOffset();
		Object pixels;
		try {
			ImageReader reader = new ImageReader(fi[0]);
			InputStream is = createInputStream(fi[0]);
			if (is == null)
				return;
			IJ.resetEscape();
			double means[] = new double[fi[0].nImages];
			for (int i = 1; i <= fi[0].nImages; i++) {
				IJ.showStatus("Reading: " + i + "/" + fi[0].nImages);
				if (IJ.escapePressed()) {
					IJ.beep();
					IJ.showProgress(1.0);
					return;
				}
				pixels = reader.readPixels(is, skip);
				if (pixels == null)
					break;
				skip = fi[0].gapBetweenImages;
				double mean = 0;
				byte[] pixels2 = (byte[]) pixels;
				for (int j = 0; j < pixels2.length; j++) {
					mean += pixels2[j] & 0xFF;
				}
				mean /= pixels2.length;
				means[i - 1] = mean;
				trace.addPoint(i, mean);
				IJ.showProgress(i, fi[0].nImages);
			}
			is.close();
			double smoothed[] = smooth(means);
			for (int i = 0; i < smoothed.length; i++) {
				filtered.addPoint(i + 1, smoothed[i]);
			}

		} catch (Exception e) {
			IJ.log("" + e);
		} catch (OutOfMemoryError e) {
			IJ.outOfMemory(fi[0].fileName);
		}
	}

	public void run(String arg) {

		// ***************** Get file ***************************

		setFile();

		setupChart();

		// Create a frame.
		frame = new JFrame("Z Axis Profile");
		Container content = frame.getContentPane();
		content.setLayout(new BoxLayout(content, BoxLayout.Y_AXIS));
		LayoutFactory factory = LayoutFactory.getInstance();
		ChartPanel chartpanel = new ChartPanel(chart);

		frame.setJMenuBar(factory.createChartMenuBar(chartpanel, false));

		content.add(chartpanel);
		content.addPropertyChangeListener(chartpanel);
		// content.add(new ControlPanel());

		// frame.add(zprofile.chart);
		frame.setVisible(true);
		// add the chart to the frame:
		frame.setSize(647, 400);
		// Enable the termination button [cross on the upper right edge]:
		JButton zoomAllButton = new JButton("Zoom All");
		zoomAllButton.addActionListener(new ZoomAllAdapter(chart));

		frame.addWindowListener(new WindowAdapter() {
			@Override
			public void windowClosing(WindowEvent e) {
				// System.exit(0);
				chart.destroy();

			}
		});

		showChart();

		plotData();

	}

	public void setFile() {
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
			this.fi = td.getTiffInfo();

			stackSize = fi[0].nImages;

		} catch (IOException e) {
			return;
		}
		if (stackSize < 2) {
			IJ.error("Need a stack");
			return;
		}
	}

	public void setupChart() {
		// ************** Create the chart ********************
		chart = new ZoomableChartWithPointSelection(this.selector);
		// Obtain the basic default axes:
		IAxis axisX = chart.getAxisX();
		IAxis axisY = chart.getAxisY();

		// Feature: Grids:
		chart.setGridColor(Color.LIGHT_GRAY);
		axisX.setPaintGrid(true);
		axisY.setPaintGrid(true);
		// Create an ITrace:
		trace = new Trace2DLtd(fi[0].nImages);
		filtered = new Trace2DSimple();
		trace.setColor(Color.GRAY);
		trace.setName("Mean intensity");
		filtered.setName("Smoothed");
		filtered.setColor(Color.blue);

		// Feature: Axis title font.
		Font titleFont = UIManager.getDefaults().getFont("Label.font")
				.deriveFont(14f).deriveFont(Font.BOLD);
		IAxis.AxisTitle axisTitle = axisY.getAxisTitle();
		axisTitle.setTitleFont(titleFont);

		// Feature: axis title.
		axisTitle.setTitle("Z Axis Profile for " + name);

		axisTitle = axisX.getAxisTitle();
		axisTitle.setTitle("slice");
		axisTitle.setTitleFont(titleFont);

		// Feature: range policy for axis.
		axisX.setRangePolicy(new RangePolicyFixedViewport(new Range(0,
				fi[0].nImages)));

		// Feature: turn on tool tips (recommended for use in static mode only):
		chart.setToolTipType(Chart2D.ToolTipType.VALUE_SNAP_TO_TRACEPOINTS);

		// Feature: turn on highlighting: Two steps enable it on the chart and
		// set a highlighter for the trace:
		Set<IPointPainter<?>> highlighters = filtered
				.getPointHighlighters();
		highlighters.clear();
		filtered.addPointHighlighter(new PointPainterDisc(12));
		chart.enablePointHighlighting(true);
		// Add all points, as it is static:
		chart.addTrace(trace);
		chart.addTrace(filtered);

	}

	public void showChart() {
		frame.setVisible(true);
	}

	public double[] smooth(double[] data) {
		// savitzky-golay coefficients
		// int[] sgcoefficients= new int[] {-21,14,39,54,59,54,39,14,-21};
		int[] sgcoefficients = new int[] { -11, 0, 9, 16, 21, 24, 25, 24, 21,
				16, 9, 0, -11 };
		double sgnorm = sum(sgcoefficients);
		int width = (int) Math.floor(sgcoefficients.length / 2);
		double smoothed[] = new double[data.length];
		for (int j = 0; j < data.length; j++) {
			smoothed[j] = 0;
			for (int k = 0; k < sgcoefficients.length; k++) {
				smoothed[j] += sgcoefficients[k]
						* data[Math.min(data.length - 1, Math.max(0, j
								- width + k))] / sgnorm;
			}
		}
		return smoothed;
	}

}
