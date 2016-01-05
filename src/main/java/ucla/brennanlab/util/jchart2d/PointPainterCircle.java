package ucla.brennanlab.util.jchart2d;

import info.monitorenter.gui.chart.ITracePoint2D;
import info.monitorenter.gui.chart.pointpainters.PointPainterDisc;

import java.awt.Graphics;

public class PointPainterCircle extends PointPainterDisc {

	/** Generated <code>serialVersionUID</code>. */
	private static final long serialVersionUID = -6317473632026920774L;

	/** The diameter of the discs to paint. */
	private int m_discSize;

	/**
	 * Cached m_discSize divided by two to save division for each point to
	 * render.
	 */
	private int m_halfDiscSize;

	/**
	 * Creates an instance with a default disc size of 4.
	 * <p>
	 */
	public PointPainterCircle() {
		this.setDiscSize(4);
	}

	/**
	 * Creates an instance with the given disc diameter.
	 * 
	 * @param diameter
	 *            the disc size in pixel to use.
	 */
	public PointPainterCircle(final int diameter) {
		this.setDiscSize(diameter);
	}

	/**
	 * Returns the diameter of the discs to paint in pixel.
	 * <p>
	 * 
	 * @return the diameter of the discs to paint in pixel.
	 */
	@Override
	public int getDiscSize() {
		return this.m_discSize;
	}

	/**
	 * @see info.monitorenter.gui.chart.IPointPainter#paintPoint(int, int, int,
	 *      int, java.awt.Graphics, info.monitorenter.gui.chart.ITracePoint2D)
	 */
	@Override
	public void paintPoint(final int absoluteX, final int absoluteY,
			final int nextX, final int nextY, final Graphics g,
			final ITracePoint2D original) {
		g.drawOval(absoluteX - this.m_halfDiscSize, absoluteY
				- this.m_halfDiscSize, this.m_discSize, this.m_discSize);
		g.fillOval(absoluteX - this.m_halfDiscSize, absoluteY
				- this.m_halfDiscSize, this.m_discSize, this.m_discSize);
	}

	/**
	 * Sets the diameter of the discs to paint in pixel.
	 * <p>
	 * 
	 * @param discSize
	 *            the diameter of the discs to paint in pixel.
	 */
	@Override
	public void setDiscSize(final int discSize) {
		this.m_discSize = discSize;
		this.m_halfDiscSize = this.m_discSize / 2;
	}

}
