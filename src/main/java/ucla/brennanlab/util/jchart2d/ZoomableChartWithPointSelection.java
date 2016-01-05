package ucla.brennanlab.util.jchart2d;

import ij.IJ;
import info.monitorenter.gui.chart.ITrace2D;
import info.monitorenter.gui.chart.ITracePoint2D;
import info.monitorenter.gui.chart.ZoomableChart;
import info.monitorenter.gui.chart.ITrace2D.DistancePoint;
import info.monitorenter.gui.chart.traces.Trace2DLtdSorted;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.util.Set;

public class ZoomableChartWithPointSelection extends ZoomableChart {

	public class ZoomAllAdapter implements ActionListener {
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

	/**
	 * 
	 */

	private double xRadius = 12; // reject new clicks if they aren't within 100
	// of previous clicks
	private static final long serialVersionUID = 3381172945559958315L;
	private boolean useSelector = false;

	// public Set<ITracePoint2D> selectedPts; // list of set points
	public ITrace2D selectedTrace;

	public ZoomableChartWithPointSelection() {
		this(false);
	}

	public ZoomableChartWithPointSelection(boolean b) {
		this.useSelector = b;
		if (!useSelector)
			return;
		new PointPainterCircle(12);
		selectedTrace = new Trace2DLtdSorted();
		selectedTrace.setTracePainter(new TracePainterCircle(12));
		selectedTrace.setColor(Color.RED);
		selectedTrace.setVisible(false);
		selectedTrace.setName("Selected points");
		selectedTrace.removeAllPoints();
		this.addTrace(selectedTrace);

	}

	@Override
	public ITracePoint2D getNearestPointManhattan(final int mouseEventX,
			final int mouseEventY) {
		ITracePoint2D result = null;
		/*
		 * Normalize pixel values:
		 */
		double scaledX = 0;
		double scaledY = 0;
		double rangeX = this.getXChartEnd() - this.getXChartStart();
		if (rangeX != 0) {
			scaledX = ((double) mouseEventX - this.getXChartStart()) / rangeX;
		}
		double rangeY = this.getYChartStart() - this.getYChartEnd();
		if (rangeY != 0) {
			scaledY = 1.0 - ((double) mouseEventY - this.getYChartEnd())
					/ rangeY;
		}

		/*
		 * TODO: Maybe cache this call because it searches all axes and evicts
		 * duplicates of their assigned traces (subject to profiling).
		 */
		Set<ITrace2D> traces = this.getTraces();
		DistancePoint distanceBean;
		DistancePoint winner = null;
		for (ITrace2D trace : traces) {
			if (!trace.isEmpty() && trace.isVisible()) {
				distanceBean = trace.getNearestPointManhattan(scaledX, scaledY);
				if (winner == null) {
					winner = distanceBean;
				} else {
					if (distanceBean.getDistance() < winner
							.getDistance()) {
						winner = distanceBean;
					}
				}
			}
		}
		if (winner != null) {
			result = winner.getPoint();
		}
		return result;
	}

	public ITracePoint2D getNearestPointManhattan(final int mouseEventX,
			final int mouseEventY, ITrace2D trace) {
		ITracePoint2D result = null;
		/*
		 * Normalize pixel values:
		 */
		double scaledX = 0;
		double scaledY = 0;
		double rangeX = this.getXChartEnd() - this.getXChartStart();
		if (rangeX != 0) {
			scaledX = ((double) mouseEventX - this.getXChartStart()) / rangeX;
		}
		double rangeY = this.getYChartStart() - this.getYChartEnd();
		if (rangeY != 0) {
			scaledY = 1.0 - ((double) mouseEventY - this.getYChartEnd())
					/ rangeY;
		}

		/*
		 * TODO: Maybe cache this call because it searches all axes and evicts
		 * duplicates of their assigned traces (subject to profiling).
		 */
		// Set<ITrace2D> traces = this.getTraces();
		DistancePoint distanceBean;
		DistancePoint winner = null;
		// for (ITrace2D trace : traces) {
		if (!trace.isEmpty() && trace.isVisible()) {
			distanceBean = trace.getNearestPointManhattan(scaledX, scaledY);
			if (winner == null) {
				winner = distanceBean;
			} else {
				if (distanceBean.getDistance() < winner
						.getDistance()) {
					winner = distanceBean;
				}
			}
		}
		// }
		if (winner != null) {
			result = winner.getPoint();
		}
		return result;
	}

	public ITracePoint2D getNearestPointManhattan(MouseEvent e, ITrace2D trace) {
		return this.getNearestPointManhattan(e.getX(), e.getY(), trace);

	}

	//

	public ITracePoint2D getNearestPointX(final int mouseEventX,
			final int mouseEventY) {
		ITracePoint2D result = null;
		/*
		 * Normalize pixel values:
		 */
		double scaledX = 0;
		double scaledY = 0;
		double rangeX = this.getXChartEnd() - this.getXChartStart();
		if (rangeX != 0) {
			scaledX = ((double) mouseEventX - this.getXChartStart()) / rangeX;
		}
		double rangeY = this.getYChartStart() - this.getYChartEnd();
		if (rangeY != 0) {
			scaledY = 1.0 - ((double) mouseEventY - this.getYChartEnd())
					/ rangeY;
		}

		/*
		 * TODO: Maybe cache this call because it searches all axes and evicts
		 * duplicates of their assigned traces (subject to profiling).
		 */
		Set<ITrace2D> traces = this.getTraces();
		DistancePoint distanceBean;
		DistancePoint winner = null;
		for (ITrace2D trace : traces) {
			if (!trace.isEmpty() && trace.isVisible()) {
				distanceBean = trace.getNearestPointManhattan(scaledX, scaledY);
				if (winner == null) {
					winner = distanceBean;
				} else {
					if (distanceBean.getDistance() < winner
							.getDistance()) {
						winner = distanceBean;
					}
				}
			}
		}
		if (winner != null) {
			result = winner.getPoint();
		}
		return result;
	}

	public ITracePoint2D getNearestPointX(final int mouseEventX,
			final int mouseEventY, ITrace2D trace) {
		ITracePoint2D result = null;
		/*
		 * Normalize pixel values:
		 */
		double scaledX = 0;
		double scaledY = 0;
		double rangeX = this.getXChartEnd() - this.getXChartStart();
		if (rangeX != 0) {
			scaledX = ((double) mouseEventX - this.getXChartStart()) / rangeX;
		}
		double rangeY = this.getYChartStart() - this.getYChartEnd();
		if (rangeY != 0) {
			scaledY = 1.0 - ((double) mouseEventY - this.getYChartEnd())
					/ rangeY;
		}

		/*
		 * TODO: Maybe cache this call because it searches all axes and evicts
		 * duplicates of their assigned traces (subject to profiling).
		 */
		// Set<ITrace2D> traces = this.getTraces();
		DistancePoint distanceBean;
		DistancePoint winner = null;
		// for (ITrace2D trace : traces) {
		if (!trace.isEmpty() && trace.isVisible()) {
			distanceBean = trace.getNearestPointManhattan(scaledX, scaledY);
			if (winner == null) {
				winner = distanceBean;
			} else {
				if (distanceBean.getDistance() < winner
						.getDistance()) {
					winner = distanceBean;
				}
			}
		}
		// }
		if (winner != null) {
			result = winner.getPoint();
		}
		return result;
	}

	public ITracePoint2D getNearestPointX(MouseEvent e, ITrace2D trace) {
		return this.getNearestPointX(e.getX(), e.getY(), trace);

	}

	@Override
	public void mouseClicked(final MouseEvent e) {
		if (e.getButton() == MouseEvent.BUTTON3) {
			return;
		}
		if (!useSelector)
			return;

		ITracePoint2D nearestpt = this.getNearestPointManhattan(e);
		;
		if (selectedTrace.isEmpty()) { // add the first point
			selectedTrace.setVisible(true);
			selectedTrace.addPoint(nearestpt);
			IJ.showStatus("You have selected a single point");

		} else {
			// was not empty before the click, now have to check the distance
			// between the new point and previous points
			// selectedTrace.addPoint(nearestpt);
			ITracePoint2D nearestselected = this.getNearestPointManhattan(e,
					selectedTrace);

			if (Math.abs(nearestselected.getX() - nearestpt.getX()) > this.xRadius) {
				selectedTrace.addPoint(nearestpt);
				IJ.showStatus("You have selected " + selectedTrace.getSize()
						+ " points");

			} else {
				IJ
						.showStatus("The point you have selected is too near other points");
				selectedTrace.removePoint(nearestselected);
				if (selectedTrace.isEmpty())
					selectedTrace.setVisible(false);
				else
					IJ.showStatus("You have selected "
							+ selectedTrace.getSize() + " points");

			}
		}

		// selectionpainter.paintPoint(e.getX(), e.getY(), 0, 0,
		// this.getGraphics(), null);

	}

}
