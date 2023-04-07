ARG TAG
ARG R_VER

FROM public.ecr.aws/chobiolab/r-dev:"$TAG"

WORKDIR /dep
COPY DESCRIPTION DESCRIPTION

ENV RENV_PATHS_LIBRARY renv/library
RUN R -e "renv::install()"

