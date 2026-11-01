FROM rocker/shiny:4.4.2

RUN R -e "install.packages(c('shiny','bslib','bsicons','DT','plotly','dplyr','tidyr','readr','janitor','ggplot2','stringr','lubridate','purrr','lme4','broom.mixed'), repos='https://cloud.r-project.org', dependencies=TRUE)"

RUN rm -rf /srv/shiny-server/*
COPY shiny_app /srv/shiny-server/
COPY data /srv/data
COPY metadata /srv/metadata

EXPOSE 3838
CMD ["/usr/bin/shiny-server"]
