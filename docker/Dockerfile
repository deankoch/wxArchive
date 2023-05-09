FROM rocker/geospatial

# copy Rscripts to working directory
ARG WX_PATH=/home/wxarchive
COPY R/* ${WX_PATH}/
WORKDIR ${WX_PATH}

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


