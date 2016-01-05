package ucla.brennanlab.imagej.util;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.measure.Measurements;
import ij.plugin.filter.Analyzer;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import info.monitorenter.gui.chart.Chart2D;
import info.monitorenter.gui.chart.IAxis;
import info.monitorenter.gui.chart.IPointHighlighter;
import info.monitorenter.gui.chart.ZoomableChart;
import info.monitorenter.gui.chart.controls.LayoutFactory;
import info.monitorenter.gui.chart.pointhighlighters.PointHighlighterConfigurable;
import info.monitorenter.gui.chart.pointpainters.PointPainterDisc;
import info.monitorenter.gui.chart.rangepolicies.RangePolicyFixedViewport;
import info.monitorenter.gui.chart.traces.Trace2DLtd;
import info.monitorenter.gui.chart.traces.Trace2DSimple;
import info.monitorenter.gui.chart.views.ChartPanel;
import info.monitorenter.util.Range;

import java.awt.Color;
import java.awt.Container;
import java.awt.Font;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.Set;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.UIManager;

import ucla.brennanlab.util.jchart2d.ZoomableChartWithPointSelection;

/** Implements the Image/Stack/Plot Z-axis Profile command. */
public class Interactive_ZAxisProfiler implements PlugInFilter, Measurements,
		ActionListener {

	class ZoomAllAdapter implements ActionListener {
		/** The zoomable chart to act upon. */
		private ZoomableChart m_zoomableChart;

		/**
		 * Creates an instance that will reset zooming on the given zoomable
		 * chart upon the triggered action.
		 * <p>
		 * 
		 * @param chart
		 *            the target to reset zooming on.
		 */
		public ZoomAllAdapter(final ZoomableChart chart) {
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

	ImagePlus imp;

	@Override
	public void actionPerformed(ActionEvent e) {
		// TODO Auto-generated method stub

	}

	@SuppressWarnings("static-access")
	float[] getZAxisProfile(Roi roi, double minThreshold, double maxThreshold) {
		ImageStack stack = imp.getStack();
		int size = stack.getSize();
		float[] values = new float[size];
		Calibration cal = imp.getCalibration();
		Analyzer analyzer = new Analyzer(imp);
		int measurements = analyzer.getMeasurements();
		boolean showResults = measurements != 0 && measurements != LIMIT;
		boolean showingLabels = (measurements & LABELS) != 0
				|| (measurements & SLICE) != 0;
		measurements |= MEAN;
		if (showResults) {
			if (!analyzer.resetCounter())
				return null;
		}
		int current = imp.getCurrentSlice();
		for (int i = 1; i <= size; i++) {
			if (showingLabels)
				imp.setSlice(i);
			ImageProcessor ip = stack.getProcessor(i);
			if (minThreshold != ImageProcessor.NO_THRESHOLD)
				ip.setThreshold(minThreshold, maxThreshold,
						ImageProcessor.NO_LUT_UPDATE);
			ip.setRoi(roi);
			ImageStatistics stats = ImageStatistics.getStatistics(ip,
					measurements, cal);
			analyzer.saveResults(stats, roi);
			if (showResults)
				analyzer.displayResults();
			values[i - 1] = (float) stats.mean;
		}
		if (showingLabels)
			imp.setSlice(current);
		return values;
	}

	public void run(ImageProcessor ip) {
		if (imp.getStackSize() < 2) {
			IJ.error("ZAxisProfiler", "This command requires a stack.");
			return;
		}
		Roi roi = imp.getRoi();
		if (roi != null && roi.isLine()) {
			IJ.error("ZAxisProfiler",
					"This command does not work with line selections.");
			return;
		}
		double minThreshold = ip.getMinThreshold();
		double maxThreshold = ip.getMaxThreshold();
		float[] y = getZAxisProfile(roi, minThreshold, maxThreshold);
		if (y != null) {
			float[] x = new float[y.length];
			for (int i = 0; i < x.length; i++)
				x[i] = i + 1;
			String title;
			if (roi != null) {
				Rectangle r = imp.getRoi().getBounds();
				title = imp.getTitle() + "-" + r.x + "-" + r.y;
			} else
				title = imp.getTitle() + "-0-0";
			/*
			 * Plot plot = new Plot(title, "Slice", "Mean", x, y); double ymin =
			 * ProfilePlot.getFixedMin(); double ymax=
			 * ProfilePlot.getFixedMax(); if (!(ymin==0.0 && ymax==0.0)) {
			 * double[] a = Tools.getMinMax(x); double xmin=a[0]; double
			 * xmax=a[1]; plot.setLimits(xmin, xmax, ymin, ymax); } plot.show();
			 */
			final ZoomableChartWithPointSelection chart = new ZoomableChartWithPointSelection();
			IAxis axisX = chart.getAxisX();
			IAxis axisY = chart.getAxisY();

			// Feature: Grids:
			chart.setGridColor(Color.LIGHT_GRAY);
			axisX.setPaintGrid(true);
			axisY.setPaintGrid(true);
			// Create an ITrace:
			Trace2DLtd trace = new Trace2DLtd(y.length);
			Trace2DSimple filtered = new Trace2DSimple();
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
			axisTitle.setTitle(title);

			axisTitle = axisX.getAxisTitle();
			axisTitle.setTitle("slice");
			axisTitle.setTitleFont(titleFont);

			// Feature: range policy for axis.
			axisX.setRangePolicy(new RangePolicyFixedViewport(new Range(0,
					y.length)));

			// Feature: turn on tool tips (recommended for use in static mode
			// only):
			chart.setToolTipType(Chart2D.ToolTipType.VALUE_SNAP_TO_TRACEPOINTS);

			// Feature: turn on highlighting: Two steps enable it on the chart
			// and
			// set a highlighter for the trace:
			Set<IPointHighlighter<?>> highlighters = filtered
					.getPointHighlighters();
			highlighters.clear();
			filtered.addPointHighlighter(new PointHighlighterConfigurable(
					new PointPainterDisc(12), true));
			chart.enablePointHighlighting(true);
			// Add all points, as it is static:
			chart.addTrace(trace);
			chart.addTrace(filtered);

			// Create a frame.
			JFrame frame = new JFrame("Z Axis Profile");
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

			double[] smoothed = smooth(y);
			for (int j = 0; j < x.length; j++) {
				trace.addPoint(x[j], y[j]);
				filtered.addPoint(x[j], smoothed[j]);
			}

			frame.addWindowListener(new WindowAdapter() {
				@Override
				public void windowClosing(WindowEvent e) {
					// System.exit(0);
					chart.destroy();

				}
			});
			frame.setVisible(true);

		}
	}

	public int setup(String arg, ImagePlus imp) {
		this.imp = imp;
		return DOES_ALL + NO_CHANGES;
	}

	public double[] smooth(float[] data) {
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
