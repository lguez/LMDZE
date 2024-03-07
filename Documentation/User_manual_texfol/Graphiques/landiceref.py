#!/usr/bin/env python3

import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs

plt.figure()
ax = plt.axes(projection=ccrs.LambertAzimuthalEqualArea(central_latitude = 90))

with xr.open_dataset("../../../../../Datasets/landiceref.nc") as f:
    f.landice.where(f.landice != 0).plot(transform = ccrs.PlateCarree(),
                                         cbar_kwargs = {"pad": 0.1})

ax.set_extent((-180,180,60,90), crs = ccrs.PlateCarree())
ax.gridlines(draw_labels=True)
##ax.coastlines()
plt.subplots_adjust(left = 0, right = 1)
plt.savefig("landiceref.pdf")
##plt.show()
