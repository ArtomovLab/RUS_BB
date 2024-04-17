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
