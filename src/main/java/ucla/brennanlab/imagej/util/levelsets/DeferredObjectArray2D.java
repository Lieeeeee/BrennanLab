package ucla.brennanlab.imagej.util.levelsets;

/**
 * Tiled array data structure for objects
 */
public class DeferredObjectArray2D<T> {
    protected final int tilesize;
    protected final Object[] tiles;
    final T defaultval;
    private final int xdim, ydim;
    private final int xtiles, ytiles;

    /**
     * Creates a new instance of DeferredDoubleArray2D
     */
    public DeferredObjectArray2D(int xdim, int ydim, int tilesize, T defaultval) {
        this.tilesize = tilesize;
        this.xdim = xdim;

        this.ydim = ydim;
        int xtiles = xdim / tilesize;
        if (xdim % tilesize > 0)
            xtiles++;
        this.xtiles = xtiles;

        int ytiles = ydim / tilesize;
        if (ydim % tilesize > 0)
            ytiles++;
        this.ytiles = ytiles;
        tiles = new Object[xtiles * ytiles];
        this.defaultval = defaultval;
    }

    public final void set(final int x, final int y, final T value) {
        final Object[][] tile = (Object[][]) getTile(x, y, true);
        tile[x % tilesize][y % tilesize] = value;
    }

    @SuppressWarnings("unchecked")
    public final T get(final int x, final int y) {
        final Object[][] tile = (Object[][]) getTile(x, y, false);

        if (tile == null) {
            return defaultval;
        } else {
            return (T) (tile[x % tilesize][y % tilesize]);
        }

    }

    public final String getAsString(final int x, final int y) {
        return this.get(x, y).toString();
    }

    protected Object createTile(final int tilesize) {
        final Object[][] tile = new Object[tilesize][tilesize];

        if (defaultval != null) {
            for (int x = 0; x < tile.length; x++) {
                for (int y = 0; y < tile[0].length; y++) {

                    tile[x][y] = defaultval;

                }
            }
        }

        return tile;
    }

    /**
     * Requests the tile for the passed coordinates.
     *
     * @param x      The X index
     * @param y      The Y index
     * @param create Determines whether the tile should be created if it has not
     *               been allocated yet.
     * @return The tile - a three dimensional array of the proper data type
     */
    protected final Object getTile(final int x, final int y,
                                   final boolean create) {
        checkBounds(x, y);

        final int x_tile = x / tilesize;
        final int y_tile = y / tilesize;

        final int offset = x_tile + y_tile * xtiles;

        final Object tile = tiles[offset];
        if (tile == null && create == true) {
            tiles[offset] = createTile(tilesize);
            return tiles[offset];
        }

        return tile;
    }

    /**
     * Returns the size of the whole virtual array in X direction
     *
     * @return The size in X direction
     */
    public final int getXLength() {
        return xdim;
    }

    /**
     * Returns the size of the whole virtual array in Y direction
     *
     * @return The size in Y direction
     */
    public final int getYLength() {
        return ydim;
    }

    private final void checkBounds(final int x, final int y) {
        if (x < 0 || x > (xdim - 1) || y < 0 || y > (ydim - 1)) {
            throw new ArrayIndexOutOfBoundsException("At index : (" + x + ", "
                    + y + ")");
        }
    }

    public int getYtiles() {
        return ytiles;
    }


}