# Downloading / Preprocessing Data

## Install DuckDB

DuckDB is used for some preprocessing steps; please follow instructions [here](https://duckdb.org/docs/installation/?version=stable&environment=cli&platform=linux&download_method=package_manager
) to install for your platform.

https://github.com/hrbrmstr/cdcfluview.git

<!-- ## ILINet Data  -->

<!-- You will need to download a dataset from [ILINet Flu Portal Dashboard](https://gis.cdc.gov/grasp/fluview/fluportaldashboard.html).  -->

<!-- On the portal, click "Download Data" and select all HHS regions and the years 2010 - 2020 in the popup; click "Download Data" on the popup to confirm. -->

<!-- ![ILINet UI](ilinet-screenshot.png) -->

<!-- Place the downloaded zip file into this directory as `FluViewPhase2Data.zip`. -->


## Running preprocessing steps

Run the command `poetry run snakemake all -j1 --snakefile prepare-data.snakefile` to complete the remaining data processing steps.





https://datarepo.eng.ucsd.edu/mcauley_group/data/amazon_v2/categoryFilesSmall/Electronics.csv
