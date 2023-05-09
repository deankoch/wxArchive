FROM rocker/geospatial

# copy all Rscripts to working directory
ARG WX_PATH=/home/wxarchive
COPY Rscript/* ${WX_PATH}/

# default arguments
ENV OPERATION=list
ENV DATA_PATH=${WX_PATH}/data

# get latest version of all required CRAN packages
RUN install2.r -e -s \
devtools \
dplyr \
jsonlite \
snapKrig \
splines \
sf \
terra \
&& rm -rf /tmp/downloaded_packages

# install latest version of my weather archive package from github
RUN R -q -e "devtools::install_github('deankoch/wxArchive')"

# run tasks with Rscript with two optional arguments
WORKDIR ${WX_PATH}
ENTRYPOINT ["sh", "-c", "Rscript run_wx.R ${OPERATION} ${DATA_PATH}"]
