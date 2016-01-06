package ucla.brennanlab.imagej.csd;

import ij.IJ;
import ij.ImagePlus;
import ij.measure.ResultsTable;
import ij.plugin.MontageMaker;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import info.monitorenter.gui.chart.ITracePoint2D;
import info.monitorenter.gui.chart.controls.LayoutFactory;
import info.monitorenter.gui.chart.views.ChartPanel;
import ucla.brennanlab.imagej.util.TIFF_Substack_Opener;
import ucla.brennanlab.imagej.util.Z_Axis_Profile_Without_Opening;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Iterator;

public class CSD_multi_montager implements PlugIn, ActionListener {
    ImagePlus montage;
    JButton montageButton;
    JFrame frame;
    int skip = 5;
    int keep = 10;
    int[] arrivalTimes;
    double[] areas;
    private Z_Axis_Profile_Without_Opening zprofile;
    private JSlider montageFrameSkip;
    private JSlider montageFrameKeep;

    @Override
    public void actionPerformed(ActionEvent e) {
        // TODO Auto-generated method stub

    }

    public void displayResultsTable() {
        if (arrivalTimes.length == 0)
            return;
        ResultsTable rt = new ResultsTable();
        for (int i = 0; i < arrivalTimes.length; i++) {
            rt.incrementCounter();

            rt.addValue("Arrival", arrivalTimes[i]);
            if (i > 0) {
                rt.addValue("Interarrival", arrivalTimes[i]
                        - arrivalTimes[i - 1]);
            }
            rt.addValue("Area", areas[i]);
        }
        rt.show("CSD arrivals");

    }

    public int[] getPeakBounds(double values[], int reference) {
        // assuming reference is on the peak, brackets the peak surrounding the
        // reference
        return null;
    }

    public void makeMontage() {
        if (zprofile.chart.selectedTrace.isEmpty())
            return;
        int numCSDs = zprofile.chart.selectedTrace.getSize();

        IJ.showStatus("Creating montages of CSDs at " + numCSDs
                + " the selected peaks...");
        montageButton.setEnabled(false);
        montageButton.setName("Creating a montage of " + numCSDs + " CSDs");
        // iterate through each of the selected points, run the csd wave peak
        // isolator, make montages
        Iterator<ITracePoint2D> iter = zprofile.chart.selectedTrace.iterator();
        int j = 1;
        TIFF_Substack_Opener subop = new TIFF_Substack_Opener(
                zprofile.directory, zprofile.name); //
        int[] stackSizes = new int[numCSDs];
        MontageMaker mon = new MontageMaker();
        ImagePlus[] peakMontages = new ImagePlus[numCSDs];
        int maxwidth = 0;
        int totalheight = 0;
        arrivalTimes = new int[numCSDs];
        areas = new double[numCSDs];
        while (iter.hasNext()) {
            IJ.showStatus("Segmenting CSD " + j + " of " + numCSDs);
            ITracePoint2D point = iter.next();
            int x = (int) point.getX();
            /***
             * 1) Find the bounds for the peak from the mean trace 2) Open an
             * ImagePlus of these bounds 3) Pass the ImagePlus into peak
             * isolator to obtain ImagePlus of new stack 4) Add result to
             * montage
             */

            ImagePlus subImage = subop.returnRange(Math.max(1, x - 100),
                    Math.min(zprofile.stackSize, x + 200));
            // subImage.show();

            CSD_peak_isolator iso = new CSD_peak_isolator(subImage, this.skip,
                    this.keep);
            iso.isolatePeaks();
            ImagePlus peakwave = iso.makeBinaryWave();
            arrivalTimes[j - 1] = iso.arrivalTime + Math.max(1, x - 100) - 1;
            areas[j - 1] = iso.area;
            peakMontages[j - 1] = mon.makeMontage2(peakwave, peakwave
                            .getStackSize(), 1, 0.5, 1, peakwave.getStackSize(), 1, 5,
                    false);
            peakMontages[j - 1].setTitle("Montage of " + zprofile.name
                    + " CSD " + j);

            if (peakMontages[j - 1].getWidth() > maxwidth)
                maxwidth = peakMontages[j - 1].getWidth();
            totalheight += peakMontages[j - 1].getHeight();

            stackSizes[j - 1] = peakwave.getStackSize();

            j++;
        }

        ImageProcessor mosaicImageIP = new ShortProcessor(maxwidth, totalheight);
        int curheight = 0;
        for (int k = 0; k < peakMontages.length; k++) {
            mosaicImageIP.insert(peakMontages[k].getProcessor(), 0, curheight);
            curheight += peakMontages[k].getHeight();
        }
        displayResultsTable();
        montage = new ImagePlus("CSD montage for " + zprofile.name,
                mosaicImageIP);
        montage.show();
        montage.updateAndDraw();
        montageButton.setName("Create Montage");
        montageButton.setEnabled(true);

    }

    public void run(String arg) {
        zprofile = new Z_Axis_Profile_Without_Opening(true);
        zprofile.setFile();
        zprofile.setupChart();

        // Setup the JFrame

        frame = new JFrame("Click on the CSD peaks to select them");

        Container content = frame.getContentPane();
        content.setLayout(new BoxLayout(content, BoxLayout.Y_AXIS));
        LayoutFactory factory = LayoutFactory.getInstance();
        ChartPanel chartpanel = new ChartPanel(zprofile.chart);

        frame.setJMenuBar(factory.createChartMenuBar(chartpanel, false));

        content.add(chartpanel);
        content.addPropertyChangeListener(chartpanel);
        content.add(new ControlPanel());

        frame.setVisible(true);
        frame.setSize(800, 600);

        zprofile.plotData();

    }

    final class ControlPanel extends JPanel {

        /**
         * Control Panel for CSD_multi_montager
         */
        private static final long serialVersionUID = 471325491258211474L;

        protected ControlPanel() {
            this.setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
            createSkipSlider();
            this.add(montageFrameSkip);
            createKeepSlider();
            this.add(montageFrameKeep);
            JComponent stretch = new JPanel();
            stretch.setBackground(Color.WHITE);
            stretch.setLayout(new BoxLayout(stretch, BoxLayout.X_AXIS));
            this.add(stretch);
            createMontageButton();
            stretch.add(Box.createHorizontalGlue());
            stretch.add(montageButton);

        }

        private void createKeepSlider() {
            montageFrameKeep = new JSlider(5, 20, keep);
            montageFrameKeep.setValue(CSD_multi_montager.this.keep);
            montageFrameKeep.setBorder(BorderFactory.createTitledBorder(
                    BorderFactory.createEtchedBorder(),
                    "Keep at most n totalSlices in the montage", TitledBorder.LEFT,
                    TitledBorder.BELOW_TOP));
            montageFrameKeep.setPaintTicks(true);
            montageFrameKeep.setPaintLabels(true);
            montageFrameKeep.setSnapToTicks(true);
            montageFrameKeep.setMajorTickSpacing(5);
            montageFrameKeep.setMinorTickSpacing(1);
            montageFrameKeep.addChangeListener(new ChangeListener() {
                public void stateChanged(final ChangeEvent e) {
                    JSlider source = (JSlider) e.getSource();
                    // Only if not currently dragged...
                    if (!source.getValueIsAdjusting()) {
                        int value = source.getValue();
                        CSD_multi_montager.this.keep = value;
                    }
                }
            });
        }

        private void createMontageButton() {
            // clear Button
            montageButton = new JButton("Make Montage");
            montageButton.setBackground(Color.WHITE);
            montageButton.setBackground(Color.WHITE);
            montageButton.addActionListener(new ActionListener() {
                public void actionPerformed(final ActionEvent e) {
                    CSD_multi_montager.this.makeMontage();
                }
            });

        }

        private void createSkipSlider() {
            montageFrameSkip = new JSlider(1, 20, skip);
            montageFrameSkip.setValue(CSD_multi_montager.this.skip);
            montageFrameSkip.setBorder(BorderFactory.createTitledBorder(
                    BorderFactory.createEtchedBorder(),
                    "Montage every n totalSlices.", TitledBorder.LEFT,
                    TitledBorder.BELOW_TOP));
            montageFrameSkip.setPaintTicks(true);
            montageFrameSkip.setPaintLabels(true);
            montageFrameSkip.setSnapToTicks(true);
            montageFrameSkip.setMajorTickSpacing(4);
            montageFrameSkip.setMinorTickSpacing(1);
            montageFrameSkip.addChangeListener(new ChangeListener() {
                public void stateChanged(final ChangeEvent e) {
                    JSlider source = (JSlider) e.getSource();
                    // Only if not currently dragged...
                    if (!source.getValueIsAdjusting()) {
                        int value = source.getValue();
                        CSD_multi_montager.this.skip = value;
                    }
                }
            });
        }
    }

}
