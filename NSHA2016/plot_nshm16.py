from matplotlib.mlab import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset as NetCDFFile
from numpy import arange, mean, percentile, array, unique, where, mgrid, ogrid,loadtxt
from matplotlib.colors import LightSource



def cpt2colormap(fileName, ncolours, **kwargs):
    import colorsys
    from numpy import array, interp, linspace
    from pylab import matplotlib

    # get kwargs
    rev = False
    for key in ['rev']:
        if key in kwargs:
            # set fault type
            if key == 'rev':
                rev = kwargs[key]

    try:
        f = open(fileName)
    except:
        print "file ", fileName, "not found"
        return None

    lines = f.readlines()
    f.close()

    x = []
    r = []
    g = []
    b = []
    colorModel = "RGB"
    for l in lines:
        ls = l.split()
        if l[0] == "#":
            if ls[-1] == "HSV":
                colorModel = "HSV"
                continue
            else:
                continue
        if ls[0] == "B" or ls[0] == "F" or ls[0] == "N":
            pass
        else:
            x.append(float(ls[0]))
            r.append(float(ls[1]))
            g.append(float(ls[2]))
            b.append(float(ls[3]))
            xtemp = float(ls[4])
            rtemp = float(ls[5])
            gtemp = float(ls[6])
            btemp = float(ls[7])

    x.append(xtemp)
    r.append(rtemp)
    g.append(gtemp)
    b.append(btemp)

    #    nTable = len(r)
    x = array(x)
    r = array(r)
    g = array(g)
    b = array(b)
    if colorModel == "HSV":
        for i in range(r.shape[0]):
            rr, gg, bb = colorsys.hsv_to_rgb(r[i] / 360., g[i], b[i])
            r[i] = rr;
            g[i] = gg;
            b[i] = bb
    if colorModel == "HSV":
        for i in range(r.shape[0]):
            rr, gg, bb = colorsys.hsv_to_rgb(r[i] / 360., g[i], b[i])
            r[i] = rr;
            g[i] = gg;
            b[i] = bb
    if colorModel == "RGB":
        r = r / 255.
        g = g / 255.
        b = b / 255.

    # reverse order
    if rev == True:
        r = r[::-1]
        g = g[::-1]
        b = b[::-1]

    # interpolate to ncolours
    xx = linspace(x[0], x[-1], ncolours)
    r = interp(xx, x, r)
    g = interp(xx, x, g)
    b = interp(xx, x, b)
    x = xx

    xNorm = (x - x[0]) / (x[-1] - x[0])

    red = []
    blue = []
    green = []
    for i in range(len(x)):
        red.append([xNorm[i], r[i], r[i]])
        green.append([xNorm[i], g[i], g[i]])
        blue.append([xNorm[i], b[i], b[i]])

    colorDict = {"red": red, "green": green, "blue": blue}

    return matplotlib.colors.LinearSegmentedColormap('my_colormap', colorDict, ncolours), xx


def get_map_polygons(m):
    polygons = []
    for polygon in m.landpolygons:
        coords = polygon.get_coords()
        tuplelist = []
        for c in coords:
            tuplelist.append((c[0], c[1]))
        polygons.append(tuplelist)

    return polygons

# for masking multiple non-conected polys
def mask_outside_polygons(polys, facecolor, plt):
    """
    Plots a mask on the specified axis ("ax", defaults to plt.gca()) such that
    all areas outside of the polygon specified by "poly_verts" are masked.
    "poly_verts" must be a list of tuples of the verticies in the polygon in
    counter-clockwise order.
    Returns the matplotlib.patches.PathPatch instance plotted on the figure.
    """
    import matplotlib.patches as mpatches
    import matplotlib.path as mpath

    # if ax is None:
    ax = plt.gca()

    # Get current plot limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    # Verticies of the plot boundaries in clockwise order
    bound_verts = [(xlim[0], ylim[0]), (xlim[0], ylim[1]),
                   (xlim[1], ylim[1]), (xlim[1], ylim[0]),
                   (xlim[0], ylim[0])]

    # A series of codes (1 and 2) to tell matplotlib whether to draw a line or
    # move the "pen" (So that there's no connecting line)
    bound_codes = [mpath.Path.MOVETO] + (len(bound_verts) - 1) * [mpath.Path.LINETO]
    verts = bound_verts

    # loop thru multiple polys
    codes = bound_codes
    for poly_verts in polys:
        verts += poly_verts[::-1]
        codes += [mpath.Path.MOVETO] + (len(poly_verts) - 1) * [mpath.Path.LINETO]

    # Plot the masking patch
    path = mpath.Path(verts, codes)
    # path = mpath.Path(vstack((bound_verts, poly_verts)), vstack((bound_codes, poly_codes)))
    patch = mpatches.PathPatch(path, facecolor=facecolor, edgecolor='none')
    patch = ax.add_patch(patch)

    # Reset the plot limits to their original extents
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    return patch

##########################################################################################
# parse grid
##########################################################################################

A = loadtxt('./NSHM16-oq29-nci/output/hazard_map-mean_3.csv',skiprows=2,delimiter=',')

# set data fields
xmllons = A[:,0]
xmllats = A[:,1]
mmi = A[:,2]

##########################################################################################
# set up map
##########################################################################################
# get map extents
urcrnrlat = 0.0
llcrnrlat = -12.
urcrnrlon = 157.
llcrnrlon = 140.5
lonspan = urcrnrlon - llcrnrlon

lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(18,10))
plt.tick_params(labelsize=16)
ax = fig.add_subplot(111)

m = Basemap(projection='merc',\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='h',area_thresh=100)

# draw coastlines, state and country boundaries, edge of map.
m.drawcoastlines()
m.drawstates()
m.drawcountries()
#m.drawmapboundary(fill_color='blue', zorder=50)
#m.fillcontinents(color='coral',lake_color='aqua')
if lonspan <= 4.0:
    tickspace = 0.5
    scale_len = 100
elif lonspan > 4.0:
    tickspace = 1.0
    scale_len = 200
m.drawparallels(arange(-90.,90.,tickspace), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.99', linewidth=0.0)
m.drawmeridians(arange(0.,360.,tickspace), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.99', linewidth=0.0)

m.drawmapscale(llcrnrlon+1, llcrnrlat+0.4, llcrnrlon+1, llcrnrlat+0.4, scale_len, \
               fontsize = 14, barstyle='fancy', zorder=100)


##########################################################################################
# read topo
##########################################################################################

mdiv = 500.

print 'Reading topo file...'
netcdffile = '../data/topo/PNG_topo.nc'
nc = NetCDFFile(netcdffile)

zscale =2. #gray
if lonspan <= 4.0:
    zscale = 30. #colour
else:
    zscale = 30. #colour

'''
# if using SRTM 30m
data = nc.variables['z'][:] / zscale
lons = nc.variables['x'][:]
lats = nc.variables['y'][:]
'''
# if using GEBCO 30 arcsec
data = nc.variables['elevation'][:] / zscale
lons = nc.variables['lon'][:]
lats = nc.variables['lat'][:]

# transform to metres
nx = int((m.xmax-m.xmin)/mdiv)+1
ny = int((m.ymax-m.ymin)/mdiv)+1
topodat = m.transform_scalar(data,lons,lats,nx,ny)

##########################################################################################
# plot intensity grid from shakemap
##########################################################################################

print 'Resampling data...'
N = 500j
extent = (llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat)

xs,ys = mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:N]
resampled = griddata(xmllons, xmllats, mmi, xs, ys, interp='linear')

# get 1D lats and lons for map transform
lons = ogrid[extent[0]:extent[1]:N]
lats = ogrid[extent[2]:extent[3]:N]

# mmidat = m.transform_scalar(resampled,lons,lats,nx,ny)
mmidat = m.transform_scalar(resampled.T, lons, lats, nx, ny)

print 'Getting colormap...'
# get colormap
cptfile = '../data/topo/haz.cpt'
cmap, zvals = cpt2colormap(cptfile, 256)


# make shading
print 'Making map...'
ls = LightSource(azdeg = 180, altdeg = 0)
norm = mpl.colors.Normalize(vmin=0.0, vmax=1.4)#myb
rgb = ls.shade_rgb(cmap(norm(mmidat)), topodat, blend_mode='hsv', vert_exag=1)
im = m.imshow(rgb)

##########################################################################################
# get land & lake polygons for masking
##########################################################################################
polys = get_map_polygons(m)


mask_outside_polygons(polys, 'lightskyblue', plt)

# get lake ploygons
polygons = []
for polygon in m.lakepolygons:
    poly = polygon.get_coords()
    plt.fill(poly[:, 0], poly[:, 1], 'lightskyblue')
    polygons.append(poly)

##########################################################################################

# set colourbar
plt.gcf().subplots_adjust(bottom=0.13)
cax = fig.add_axes([0.34,0.05,0.33,0.03]) # setup colorbar axes.
cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')

plt.show()
