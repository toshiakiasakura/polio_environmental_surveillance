FROM jupyter/datascience-notebook:2023-01-09

WORKDIR /workdir

ENV JULIA_PROJECT=/workdir
ENV JULIA_DEPOT_PATH=":/opt/julia"

RUN rm -rf /opt/julia/registries && \
    mkdir /opt/julia/registries && \
    git clone https://github.com/JuliaRegistries/General.git /opt/julia/registries

RUN pip install jupytext==1.14.7
