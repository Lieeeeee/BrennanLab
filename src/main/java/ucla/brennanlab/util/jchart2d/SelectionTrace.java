package ucla.brennanlab.util.jchart2d;

import info.monitorenter.gui.chart.traces.Trace2DLtdSorted;

public class SelectionTrace extends Trace2DLtdSorted {
	/**
	 * 
	 */
	private static final long serialVersionUID = -4719868771114763860L;
	boolean tooltip = true;

	public SelectionTrace() {
		super(100);
	}

	public SelectionTrace(int size) {
		super(size);
	}
}