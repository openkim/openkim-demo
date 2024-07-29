# base image
FROM ghcr.io/openkim/developer-platform@sha256:d30f0c5d0164196d0df89b8daf1189da39218696c0ff736d5579d387159e0fb6

USER root
RUN pip install --no-cache-dir notebook jupyterlab nglview

# Make sure the contents of our repo are in ${HOME}
COPY . /home/openkim/openkim-demo/
RUN chown -R openkim /home/openkim
USER openkim
ENV SHELL=/bin/bash
ENV LD_LIBRARY_PATH=:/usr/local/lib
