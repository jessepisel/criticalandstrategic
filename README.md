# Wyoming Critical and Strategic Minerals

This repository is home to the data processing scripts for processing [NURE](https://mrdata.usgs.gov/nure/sediment/) and [NGDBR](https://mrdata.usgs.gov/ngdb/rock/) geochemical datasets as described in [WSGS OFR 2019-2](http://www.wsgs.wyo.gov/products/wsgs-2019-ofr-02.pdf)

All scripts were written in Jupyter notebooks, and all imports are stated at the top of each notebook. Notebooks are listed in order of use from 01-04. The only calculation outside of Jupyter notebooks is calculating the *Gi** value for the rock sample dataset. This was done in ArcMap, and exported to `xls` workbooks.

## Notebooks are as follows
01 Imputation.ipynb - takes rock chip and sediment sample geochemistry data in `CSV` format and imputes missing values. It then takes the data and selects the highest value for each element for each drainage basin and creates a shapefile. 

02 Imputation validation.ipynb - compares imputed values from the first notebook to actual measured values and creates visualizations

03 Sediment hotspots.ipynb - calculates the *Gi** value for sediment samples using a dendritic spatial weights matrix

04 Hotspots to classes.ipynb - combines the rock and sediment hotspot data for classification
