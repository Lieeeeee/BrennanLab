package ucla.brennanlab.imagej.util.levelsets;


/**
 * Voxel element for the Fast Marching trial set.
 */
public class FastMarchingBandElement extends BandElement
{
   // The distance value
   private double distance = 0;
   
   /**
    * Creates a new instance of FastMarchingBandElement
    * @param x The X coordinate
    * @param y The Y coordinate
    * @param value The X coordinate
    */
   public FastMarchingBandElement(int x, int y, int value)
   {
      super(x, y , value);
   }

   /**
    * Returns the distance value
    * @return The distance value
    */
   public double getDistance()
   {
      return distance;
   }

   /**
    * Sets the distance value
    * @param distance The distance value to be set
    */
   public void setDistance(double distance)
   {
      this.distance = distance;
   }  
}