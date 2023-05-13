FROM rocker/geospatial

# working directory in container
ENV WX_PATH=/home/wxarchive

# mount your local storage volume to this location
ENV WX_DATA_DIR=/home/wxarchive/data

# area of interest polygon
ARG AOI_NAME=aoi.geojson

# get CRAN packages (skip any installed already)
RUN install2.r -e -s \
devtools \
dplyr \
jsonlite \
snapKrig \
splines \
sf \
terra \
&& rm -rf /tmp/downloaded_packages

# copy the AOI file and install wxArchive
WORKDIR ${WX_PATH}
COPY ${AOI_NAME} ${AOI_NAME}
RUN R -q -e "devtools::install_github('deankoch/wxArchive', force=TRUE)"

# copy Rscript(s)
COPY Rscript/* ${WX_PATH}/

# set default operation and run Rscript
ENV WX_OPERATION=list
ENTRYPOINT ["sh", "-c", "Rscript ${WX_PATH}/run_wx.R ${WX_OPERATION} ${WX_DATA_DIR}"]
