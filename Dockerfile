FROM jupyter/datascience-notebook:2023-01-09

WORKDIR /workdir

#RUN R -e "library(remotes); install_github('stewid/SimInf@v9.5.0')"
RUN pip install jupyterlab-vim

ENV JULIA_PROJECT=/workdir
RUN julia -e 'using Pkg; Pkg.activate(".")'
COPY ./Manifest.toml Manifest.toml
COPY ./Project.toml Project.toml
RUN julia -e 'using Pkg; Pkg.instantiate()'
