FROM rocker/geospatial

# working directory in container
ENV WX_PATH=/home/wxarchive

# mount your local storage volume to this location
ENV WX_DATA_DIR=/home/wxarchive/data

# area of interest polygon
ARG AOI_NAME=aoi.geojson
ARG EXP_NAME=aoi_export.geojson

# fitted temporal model files
ARG MODEL_NAME=model.zip

# get CRAN packages (skip any installed already)
RUN install2.r -e -s \
devtools \
dplyr \
jsonlite \
snapKrig \
splines \
sf \
s2 \
terra \
&& rm -rf /tmp/downloaded_packages

# copy the AOI and model files, and install wxArchive
WORKDIR ${WX_PATH}
COPY ${AOI_NAME} ${AOI_NAME}
COPY ${EXP_NAME} ${EXP_NAME}
COPY ${MODEL_NAME} ${MODEL_NAME}
RUN R -q -e "devtools::install_github('deankoch/wxArchive', force=TRUE)"

# copy Rscript(s)
COPY Rscript/* ${WX_PATH}/

# WX_START and WX_END can be strings like "2023-01-01"
ENV WX_OPERATION=list
ENV WX_START='auto'
ENV WX_END='auto'

# run the Rscript that was copied above
ENTRYPOINT ["sh", "-c", "Rscript ${WX_PATH}/run_wx.R ${WX_DATA_DIR} ${WX_OPERATION} ${WX_START} ${WX_END}"]
