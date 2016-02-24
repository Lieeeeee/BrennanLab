package org.brennanlab.util.jchart2d;

import info.monitorenter.gui.chart.ITracePoint2D;
import info.monitorenter.gui.chart.pointpainters.PointPainterDisc;
import info.monitorenter.gui.chart.traces.painters.ATracePainter;

import java.awt.*;

/**
 * Renders traces by painting a disc (hollow circle) with choosable diameter for
 * each {@link info.monitorenter.gui.chart.TracePoint2D} to show.
 * <p>
 *
 * @author <a href="mailto:Achim.Westermann@gmx.de">Achim Westermann </a>
 * @version $Revision: 1.17 $
 */
public class TracePainterCircle extends ATracePainter {

    /**
     * Generated <code>serialVersionUID</code>.
     */
    private static final long serialVersionUID = 8919406018882664083L;

    /**
     * The implementation for rendering the point as a disc.
     */
    private PointPainterDisc m_pointPainter;

    /**
     * Creates an instance with a default disc size of 4.
     * <p>
     */
    public TracePainterCircle() {
        this.m_pointPainter = new PointPainterCircle(4);
    }

    /**
     * Creates an instance with the given disc size.
     *
     * @param discSize the disc size in pixel to use.
     */
    public TracePainterCircle(final int discSize) {
        this.m_pointPainter = new PointPainterCircle(discSize);
    }

    /**
     * @see info.monitorenter.gui.chart.ITracePainter#endPaintIteration(java.awt.Graphics)
     */
    @Override
    public void endPaintIteration(final Graphics g2d) {
        if (g2d != null) {
            this.m_pointPainter.paintPoint(this.getPreviousX(), this
                    .getPreviousY(), 0, 0, g2d, this.getPreviousPoint());
        }
        this.m_pointPainter.endPaintIteration(g2d);
    }

    /**
     * Returns the diameter of the discs to paint in pixel.
     * <p>
     *
     * @return the diameter of the discs to paint in pixel.
     */
    public int getDiscSize() {
        return this.m_pointPainter.getDiscSize();
    }

    /**
     * Sets the diameter of the discs to paint in pixel.
     * <p>
     *
     * @param discSize the diameter of the discs to paint in pixel.
     */
    public void setDiscSize(final int discSize) {
        this.m_pointPainter.setDiscSize(discSize);
    }

    /**
     * @see info.monitorenter.gui.chart.traces.painters.ATracePainter#paintPoint(int,
     * int, int, int, java.awt.Graphics,
     * info.monitorenter.gui.chart.ITracePoint2D)
     */
    @Override
    public void paintPoint(final int absoluteX, final int absoluteY,
                           final int nextX, final int nextY, final Graphics g,
                           final ITracePoint2D original) {
        super.paintPoint(absoluteX, absoluteY, nextX, nextY, g, original);
        this.m_pointPainter.paintPoint(absoluteX, absoluteY, nextX, nextY, g,
                original);
    }

    /**
     * @see info.monitorenter.gui.chart.traces.painters.ATracePainter#startPaintIteration(java.awt.Graphics)
     */
    @Override
    public void startPaintIteration(final Graphics g2d) {
        this.m_pointPainter.startPaintIteration(g2d);
    }
}
