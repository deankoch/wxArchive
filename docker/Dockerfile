FROM rocker/geospatial

# working directory in container
ENV WX_PATH=/home/wxarchive
WORKDIR ${WX_PATH}

# mount your local storage volume to this location
ENV WX_DATA_DIR=/home/wxarchive/data

# source tarball for wxArchive
ENV PKG_NAME=wxArchive_0.0.0.9000.tar.gz
COPY ${PKG_NAME} ${PKG_NAME}

# area of interest polygon
ARG AOI_NAME=aoi.geojson
COPY ${AOI_NAME} ${AOI_NAME}

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

# install wxArchive
RUN R -q -e "install.packages('${WX_PATH}/${PKG_NAME}', repo=NULL, type='source')"

# copy Rscript(s)
COPY Rscript/* ${WX_PATH}/

# set default operation and run Rscript
ENV WX_OPERATION=list
ENTRYPOINT ["sh", "-c", "Rscript ${WX_PATH}/run_wx.R ${WX_OPERATION} ${WX_DATA_DIR}"]
