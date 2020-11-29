# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 21:41:44 2020

@author: patry
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 14:49:02 2020

@author: patry
"""
import numpy as np
from astropy.table import Table, join
from matplotlib import pyplot as plt
import astropy.io.fits as pyfits
import os
from astropy.stats import binom_conf_interval


datapath = "C:/Users/patry/OneDrive/Documents/Python Scripts"
# I don't have your matched table, so I'll have to do it myself...

gkgz_table = Table.read(os.path.join(datapath, "gkgz-debiased-final.fits"))
env_table = Table.read(os.path.join(datapath, "EnvironmentMeasuresv05.fits"))
mass_table = Table.read(os.path.join(datapath, "StellarMassesLambdarv20.fits"))

print(len(gkgz_table), len(env_table), len(mass_table), )

gkgz_mass_table = join(gkgz_table, mass_table, keys='CATAID', metadata_conflicts='silent')
print(len(gkgz_mass_table))

d = join(gkgz_mass_table, env_table, keys='CATAID', metadata_conflicts='silent')
print(len(d))

duplicate_columns = [c for c in d.colnames if c.endswith('_1')]
for col1 in duplicate_columns:
    col2 = col1.replace('_1', '_2')
    col0 = col1.replace('_1', '')
    match = (d[col1] - d[col2] < 1e-3).all()
    if match:
        d.rename_column(col1, col0)
        d.remove_column(col2)
        print(f'Matching columns combined in {col0}')
    else:
        print(f'Duplicate columns {col1}, {col2} do not match')
        
# setting negative spiral frac values to 0
# (this includes both galaxies without features and edge-on galaxies)
d['spiral_spiral_deb_frac'][d['spiral_spiral_deb_frac'] < 0] = 0
d['spiral'] = d['spiral_spiral_deb_frac'] > 0.5
# calculate disk galaxies as those with spiral arms or endge-on
d['disk'] = d['spiral_spiral_deb_frac'] > 0.5
d['disk'] |= d['edgeon_yes_deb_frac'] > 0.5

#plt.figure()
#counts, bins, p = plt.hist(d['features_features_deb_frac'])
#plt.hist(d['spiral_spiral_deb_frac'], bins=bins, histtype='step', lw=2);
#plt.figure()
#plt.hist(d['Z_TONRY'], bins=50);

# One or two galaxies have 'bad' absmag values
ok = d['absmag_r'] > -90

# Restrict to only reliably measured surface densities
ok &= d['SurfaceDensityFlag'] == 0

full_table = d[ok]

z_min = 0.03
z_lim = 0.085
absmag_lim = -20.17
#absmag_lim = -18  # we could go fainter with this dataset

#z_lim = 0.1  # or further
#absmag_lim = -18.5

#z_lim = 0.15
#absmag_lim = -19.25

in_vol_limit = full_table['Z_TONRY'] > z_min
in_vol_limit &= full_table['Z_TONRY'] < z_lim
in_vol_limit &= full_table['absmag_r'] < absmag_lim


#plt.figure()
vol_table = full_table[in_vol_limit]
vol_table['logSurfaceDensity'] = np.log10(vol_table['SurfaceDensity'])
c1, b, p= plt.hist(vol_table['logSurfaceDensity'], bins = 15);

limit = vol_table['logSurfaceDensity'] > -0.5
limit &= vol_table['logSurfaceDensity'] < 0.65
vol_table = vol_table[limit]

plt.hist(vol_table['logSurfaceDensity'], bins = b,histtype='step', lw=2)


#GROUPING BY MASS

counts, bins = np.histogram(vol_table['logmstar'], bins=10, range=[9.5,11.3])
vol_table['mass_bins']= np.digitize(vol_table['logmstar'], bins)
grouped = vol_table.group_by('mass_bins')

fig, ax = plt.subplots()

plt.xlabel('log10[SurfaceDensity]')
plt.ylabel('Spiral Fraction')
plt.axis()
legend = ax.legend(loc='upper right', shadow=True, prop={'size': 6})

#environmental bins
#for all mass bins group by galaxy density
for x in range(1, 10):
    mb = grouped['mass_bins'] == x
    vol_table1 = grouped[mb] 
    lens=len(vol_table1['mass_bins'])
    
    #only consider bins with data points more than 100
    if lens > 100:
        r1 = vol_table1['logSurfaceDensity'].min()
        r2 = vol_table1['logSurfaceDensity'].max()
        count1, bins1 = np.histogram(vol_table1['logSurfaceDensity'], bins=5, range=[r1,r2])
        vol_table1['env_bins']= np.digitize(vol_table1['logSurfaceDensity'], bins1)
        grouped1 = vol_table1.group_by('env_bins')
        means = grouped1.groups.aggregate(np.mean)
        
        #creating an array to find number of data points in each environmental bin
        n=[]
        #thios doesnt work properly, n ends up having different dimensions to p 
        for y in range(1,6):
            lens1= vol_table1['env_bins'] == y
            yy = grouped1[lens1]
            yy =len(yy['env_bins'])
            n = n + [yy]
            
        #binomial errors   
        p= means['spiral_spiral_deb_frac']
        
        k = n*p
        a, b= binom_conf_interval(k,n)
        error = [a,b]
        #unsure what to set the upper and lower limits as 
        ax.errorbar(means['logSurfaceDensity'], means['spiral_spiral_deb_frac'],
                error, lolims = p/2, uplims = p/2, fmt='.-')



