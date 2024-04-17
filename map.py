%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')
from mpl_toolkits.basemap import Basemap
from itertools import chain

def draw_map(m, scale=0.2):

    m.shadedrelief(scale=scale)
    
    lats = m.drawparallels(np.linspace(-90, 90, 13))
    lons = m.drawmeridians(np.linspace(-180, 180, 13))

    lat_lines = chain(*(tup[1][0] for tup in lats.items()))
    lon_lines = chain(*(tup[1][0] for tup in lons.items()))
    all_lines = chain(lat_lines, lon_lines)
    
    for line in all_lines:
        line.set(linestyle='-', alpha=0.3, color='w')

fig = plt.figure(figsize=(40, 40))
m = Basemap(projection='lcc', resolution='l',
            lon_0=93, lat_0=67, lat_1=60, lat_2=80,
            width=0.9E7, height=0.5E7)

m.shadedrelief()
draw_map(m)
plt.show()
plt.savefig('map_RUS.png')


fig = plt.figure(figsize=(20, 20))
coord = pd.read_table('coordinates.txt',delimiter='\t')
map = Basemap(projection='cyl',llcrnrlat=-60,urcrnrlat=90,
            llcrnrlon=-168,urcrnrlon=192, resolution='l')
map.drawmapboundary(fill_color='#FFFFFF')
map.fillcontinents(color='#afe9c6ff',lake_color='#FFFFFF')
map.drawcoastlines(linewidth=0.25, linestyle='solid',antialiased=1)


cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["grey","yellow","red"])
map.scatter(coord['lon'], coord['lat'], c=coord['IBD_mean'], s=250, cmap=cmap,marker='o',zorder=10000)

#plt.legend(handles=plot.legend_elements()[0])
plt.colorbar()

plt.savefig('map_add_1000G.pdf',format = 'pdf')

#for several populations coordinates were modified
'''
G_coord[G_coord$POP2 == 'CEU',]$lon <- 2.714358
G_coord[G_coord$POP2 == 'CEU',]$lat <- 47.422519

G_coord[G_coord$POP2 == 'STU',]$lon <- 80.71378
G_coord[G_coord$POP2 == 'STU',]$lat <- 7.555494

G_coord[G_coord$POP2 == 'ITU',]$lon <- 79.11517
G_coord[G_coord$POP2 == 'ITU',]$lat <- 17.84959

G_coord[G_coord$POP2 == 'GIH',]$lon <- 71.745263
G_coord[G_coord$POP2 == 'GIH',]$lat <- 22.385001
'''
