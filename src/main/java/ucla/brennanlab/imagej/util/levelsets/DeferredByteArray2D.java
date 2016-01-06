package ucla.brennanlab.imagej.util.levelsets;

public class DeferredByteArray2D {
    protected final int tilesize;
    protected final Object[] tiles;
    private final int xdim, ydim;
    private final int xtiles, ytiles;
    byte defaultval = 0;

    public DeferredByteArray2D(final int xdim, final int ydim,
                               final int tilesize, final byte defaultval) {
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
    }

    public final int getTileSize() {
        return tilesize;
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

    /**
     * Creates a tile. This is delegated to concrete subclasses as the tile
     * needs to be of the proper data type.
     *
     * @param tilesize The tile dimension
     * @return The tile - a three dimensional array of the proper data type
     */
    protected final Object createTile(final int tilesize) {
        final byte[][] tile = new byte[tilesize][tilesize];

        if (defaultval != 0) {
            for (int x = 0; x < tile.length; x++) {
                for (int y = 0; y < tile[0].length; y++) {

                    tile[x][y] = defaultval;
                }
            }
        }

        return tile;
    }

    public final byte get(final int x, final int y) {
        final byte[][] tile = (byte[][]) getTile(x, y, false);

        if (tile == null) {
            return defaultval;
        } else {
            return tile[x % tilesize][y % tilesize];
        }

    }

    public final void set(final int x, final int y, final byte value) {
        final byte[][] tile = (byte[][]) getTile(x, y, true);
        tile[x % tilesize][y % tilesize] = value;
    }

    public int getYtiles() {
        return ytiles;
    }

}