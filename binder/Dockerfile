# base image
FROM ghcr.io/openkim/developer-platform@sha256:245afa47aa93fc9729616d82999dad7ae707a05f7a7e03f57b6a2fdc7f14bc76

USER root
RUN pip install --no-cache-dir notebook jupyterlab==4.3.6 nglview==3.1.2

# Make sure the contents of our repo are in ${HOME}
COPY . /home/openkim/openkim-demo/
RUN chown -R openkim /home/openkim
USER openkim
ENV SHELL=/bin/bash
