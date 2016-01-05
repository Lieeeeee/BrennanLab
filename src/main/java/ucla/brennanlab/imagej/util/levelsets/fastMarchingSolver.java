package ucla.brennanlab.imagej.util.levelsets;

import ij.IJ;
import java.util.PriorityQueue;

public class fastMarchingSolver {
	// computation matrix and vectors for fast marching

	private PriorityQueue<BandElement> heap = null;
	private BandElementCache elem_cache = null;

	private ImplicitShape2D initial; // Initial location at t=0;
	// Default parameters

	// Maximum distance to be traveled (pixels on a path from seed-point)
	private static final double DISTANCE_STOP = 1000; // 1;
	// Growth threshold for immediate abort
	
	private static final int ERROR_STOP = 5;
	// stop if there are a certain number of root finding errors
	
	private static final double EXTREME_GROWTH = 1000;
	// Arrival time of the last voxel added to the alive set
	private double lastFreezeTime = 0;
	// Maximum distance in pixels from seed point travelled so far
	private double max_distance = 0;

	// Holds the state of the voxels
	private DeferredByteArray2D map = null;
	// Holds the arrival timeof the voxels
	private DeferredDoubleArray2D arrival = null;
	/*
	 * Holds the distance (shortest way over already visited voxels) from seed
	 * point of the voxels
	 */
	private DeferredDoubleArray2D distances = null;
	private DeferredObjectArray2D<BandElement> elementLUT = null;

	private final static int ELEMENT_CACHE_SIZE = 1000;
	// control parameters for fast marching
	float[][] spd;
	int width, height;

	/**
	 * Constant for Far elements
	 */
	public static final byte FAR = 0;
	/**
	 * Constant for Band elements
	 */
	public static final byte BAND = 1;
	/**
	 * Constant for ALIVE elements
	 */
	public static final byte ALIVE = 2;
	
	public fastMarchingSolver(ImplicitShape2D initial){
		this.setInitial(initial);
	}

	/**
	 * Initialize fast marching using a given level set position and speed
	 * 
	 * @param initial
	 * @param spd
	 */
	public fastMarchingSolver(ImplicitShape2D initial, float[][] spd) {
		this.setInitial(initial);
		this.spd=spd;
		reset();

	}
	
	public void setSpd(float[][] spd){
		this.spd=spd;
	}
	
	
	public void reset(){
		heap = new PriorityQueue<BandElement>(1000);
		elem_cache = new BandElementCache(ELEMENT_CACHE_SIZE);
		map = new DeferredByteArray2D(initial.width, initial.height, 5, FAR);
		arrival = new DeferredDoubleArray2D(initial.width, initial.height, 5, 0);
		distances = new DeferredDoubleArray2D(initial.width, initial.height, 5,  0d);    
		elementLUT = new DeferredObjectArray2D<BandElement>(initial.width, initial.height, 5,null);
				
		// First make sure we actually got valid seed points
		if (initial == null) {
			throw new IllegalArgumentException(
					"Fast Marching needs seed points but didn't find any! Did you specify an area?");

		}
		
		for(int x=0;x<initial.width;x++){
			for(int y=0;y<initial.height;y++){
				if(Math.abs(initial.get(x, y))<=1 && initial.get(x, y)>-1){
					final Coordinate seed = new Coordinate(x,y);
					final BandElement start = elem_cache.getRecycledBandElement(seed.getX(),seed.getY(),0);
					elementLUT.set(seed.getX(), seed.getY(), start);
					map.set(seed.getX(), seed.getY(), BAND);
					heap.add(start);
				}
				if(initial.get(x,y)<1){
					final Coordinate seed = new Coordinate(x,y);
					final BandElement start = elem_cache.getRecycledBandElement(seed.getX(),seed.getY(),0);
					elementLUT.set(seed.getX(), seed.getY(), start);
					map.set(seed.getX(), seed.getY(), ALIVE);
				}
			}

			
		}
	}
	
	
	public ImplicitShape2D maskToLevelSet(DeferredByteArray2D mask){
		boolean[][] boolMask = new boolean[mask.getXLength()][mask.getYLength()];
		for(int x = 0; x<mask.getXLength();x++){
			for(int y=0;y<mask.getYLength();y++){
				boolMask[x][y] = mask.get(x, y) == ALIVE ? true : false;
			}
		}
		return new ImplicitShape2D(boolMask);
	}
	
	/**
	 * Solve forward in time to either time or to end of spatial domain
	 * @param time
	 * @return Level set representation of the interface
	 */
	public ImplicitShape2D solveAndReturnSignedDistance(double time){
		
		if(heap.isEmpty()){
			cleanup();
			return maskToLevelSet(this.map);
		}
		double t=0;
		
		while(t<time){
			final BandElement next = heap.poll();
			
			freeze(next);
			t=arrival.get(next.getX(), next.getY());
			elem_cache.recycleBandElement(next);
			if(heap.isEmpty()){
				cleanup();
				return  maskToLevelSet(this.map);
			}
		}
		
		return maskToLevelSet(this.map);
	}
	// dereference large data structure to allow garbage collection
	     private final void cleanup()
	     {
	        arrival = null;
	        elementLUT = null;
	        heap = null;
	        elem_cache = null;
	     }
	private final void freeze(final BandElement elem) {
		final int freezeX = elem.getX();
		final int freezeY = elem.getY();
		map.set(freezeX, freezeY, ALIVE);
		elementLUT.set(freezeX, freezeY, null);

		final double dist = distances.get(freezeX, freezeY);

		if (dist > DISTANCE_STOP) {
			System.out.println("Stopped - max distance exceeded");
			heap.clear();
			return;
		}/* else if(false)// (dist < max_distance )

		{
			arrival.set(freezeX, freezeY, Double.MAX_VALUE);
			IJ.log("whyyyyyyyyy");
			distances.set(freezeX, freezeY, Double.MAX_VALUE);
			System.out.println("Sorted out voxel");
		} */else {
			arrival.set(freezeX, freezeY, elem.getValue());

			if (max_distance < distances.get(freezeX, freezeY)) {
				max_distance = distances.get(freezeX, freezeY);
			}

			if (freezeX > 0)
				update(freezeX - 1, freezeY);
			if (freezeX + 1 < map.getXLength())
				update(freezeX + 1, freezeY);
			if (freezeY > 0)
				update(freezeX, freezeY - 1);
			if (freezeY + 1 < map.getYLength())
				update(freezeX, freezeY + 1);

			if (arrival.get(freezeX, freezeY) > (lastFreezeTime + EXTREME_GROWTH)) {

				IJ.log("Fast marching stopped - extreme growth");
				IJ.log("Last -> " + lastFreezeTime);
				IJ.log("Now -> " + arrival.get(freezeX, freezeY));
				heap.clear();
				elementLUT = null;
			} else {
				lastFreezeTime = arrival.get(freezeX, freezeY);
			}
			// }
		}

	}

	// Updates a voxel in the neighborhood of voxel just moved to alive set
	private final void update(final int x, final int y) {
		final byte cell_state = map.get(x, y);
		if (cell_state == ALIVE)
			return;

		final double time = calculateArrivalTime(x, y);
		final double dist = calculateDistance(x, y);

		// If this voxel is already in the trial update arrival time an dd
		// distance
		if (cell_state == BAND) {
			final BandElement elem = elementLUT.get(x, y);

			/*
			 * updated distance and arrival time is guaranteed to be <= old
			 * distance so omit a time consuming check
			 */

			heap.remove(elem);
			elem.setValue(time);
			heap.offer(elem);

			distances.set(x, y, dist);
		}
		// If this voxel is currently in the far set add it to the trial set
		else if (cell_state == FAR) {
			final BandElement elem = elem_cache.getRecycledBandElement(x, y,
					time);
			heap.offer(elem);

			map.set(x, y, BAND);
			elementLUT.set(x, y, elem);
			distances.set(x, y, dist);
		}
	}

	private final double calculateArrivalTime(final int x, final int y) {
		// Get neighbor with minimal arrival time in every spatial direction

		final double xB = (x > 0 && map.get(x - 1, y) == ALIVE) ? arrival.get(
				x - 1, y) : Double.MAX_VALUE;
		final double xF = (x + 1 < map.getXLength() && map.get(x + 1, y) == ALIVE) ? arrival
				.get(x + 1, y)
				: Double.MAX_VALUE;
		final double yB = (y > 0 && map.get(x, y - 1) == ALIVE) ? arrival.get(
				x, y - 1) : Double.MAX_VALUE;
		final double yF = (y + 1 < map.getYLength() && map.get(x, y + 1) == ALIVE) ? arrival
				.get(x, y + 1)
				: Double.MAX_VALUE;

		final double xVal = (xB < xF) ? xB : xF;
		final double yVal = (yB < yF) ? yB : yF;

		// Determine quadratic coefficient.
		int quadCoeff = 0;
		if (xVal < Double.MAX_VALUE)
			quadCoeff++;
		if (yVal < Double.MAX_VALUE)
			quadCoeff++;

		final float speed = this.spd[x][y];

		/*
		 * If only one spatial direction contributes to the quadratic
		 * coefficient, than there ist a much more efficient solution - so this
		 * is calculated and returned instead then.
		 */
		if (quadCoeff == 1) {
			if (xVal < Double.MAX_VALUE) {
				return xVal + (1 / speed);
			} else if (yVal < Double.MAX_VALUE) {
				return yVal + (1 / speed);
			}
		}

		double solution = 0;
		double linCoeff = 0;
		double abs = (-1 / (speed * speed));

		/*
		 * Calculate linear and absolute term contributions for every spatial
		 * direction if alive.
		 */
		if (xVal < Double.MAX_VALUE) {
			linCoeff -= 2 * xVal;
			abs += xVal * xVal;
		}
		if (yVal < Double.MAX_VALUE) {
			linCoeff -= 2 * yVal;
			abs += yVal * yVal;
		}

		// Discriminant of the general quadratic equation
		final double discriminant = (linCoeff * linCoeff)
				- (4 * quadCoeff * abs);

		// Two solutions exist. Calculate the bigger one.
		if (discriminant > 0) {
			final double rootDiscriminant = Math.sqrt(discriminant);
			solution = ((-linCoeff) + rootDiscriminant) / (2 * quadCoeff);

		}
		// No solution exists - read below
		else {
			/*
			 * Something went really wrong - no solution for the quadratic
			 * equation. This should NEVER happen, so it clearly indicates a
			 * problem with the speed calculation.
			 */
			IJ.log("OUCH !!! # solutions = 0 (at " + x + ", " + y + ")");
			IJ.log("quad. coefficient = " + quadCoeff);
			IJ.log("lin. coefficient = " + linCoeff);
			IJ.log("absolute term = " + abs);
			IJ.log("xVal = " + xVal);
			IJ.log("yVal = " + yVal);
			IJ.log("Speedterm = " + speed);
			IJ
					.log("xB, xF, yB, yF = " + xB + ", " + xF + ", " + yB
							+ ", " + yF);
			IJ.log("**********************************");
			return Double.NaN;
			//System.exit(0);
		}

		return solution;
	}

	// Calculate the distance to the nearest seed-point using only alive
	// waypoints
	private final double calculateDistance(final int x, final int y) {
		// Get distances of all alive neighbors
		final double xB = (x > 0 && map.get(x - 1, y) == ALIVE) ? distances
				.get(x - 1, y) : Double.MAX_VALUE;
		final double xF = (x + 1 < map.getXLength() && map.get(x + 1, y) == ALIVE) ? distances
				.get(x + 1, y)
				: Double.MAX_VALUE;
		final double yB = (y > 0 && map.get(x, y - 1) == ALIVE) ? distances
				.get(x, y - 1) : Double.MAX_VALUE;
		final double yF = (y + 1 < map.getYLength() && map.get(x, y + 1) == ALIVE) ? distances
				.get(x, y + 1)
				: Double.MAX_VALUE;

		// Find minimum of the distances
		final double xVal = (xB < xF) ? xB : xF;
		final double yVal = (yB < yF) ? yB : yF;

		final double dist = Math.min(xVal, yVal);

		// Add 1 to the smallest way of its neighbors to reach this voxel
		return (dist + 1);
	}


	public void setInitial(ImplicitShape2D initial) {
		this.initial = initial;
	}


	public ImplicitShape2D getInitial() {
		return initial;
	}


	public boolean isDone() {
		if(heap==null) return true;
		if(heap.isEmpty()) return true;
		return false;
	}

	public static int getErrorStop() {
	    return ERROR_STOP;
	}

}
