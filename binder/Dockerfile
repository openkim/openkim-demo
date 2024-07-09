# base image
FROM ghcr.io/openkim/developer-platform@sha256:3623be2f82cb7a00df735dbda00bffc5d7395862ccd41528df7f33884548da77

USER root
RUN pip install --no-cache-dir notebook jupyterlab nglview

# Make sure the contents of our repo are in ${HOME}
COPY . /home/openkim/openkim-demo/
RUN chown -R openkim /home/openkim
USER openkim
ENV SHELL=/bin/bash
ENV LD_LIBRARY_PATH=:/usr/local/lib
