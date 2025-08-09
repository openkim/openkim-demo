# base image
FROM ghcr.io/openkim/developer-platform@sha256:9f5596c4bfb60426b71e628f1a4f4798ad1b9656585b7d3464edee99bac1751e

USER root
RUN pip install --no-cache-dir notebook jupyterlab==4.3.6 nglview==3.1.2

# Make sure the contents of our repo are in ${HOME}
COPY . /home/openkim/openkim-demo/
RUN chown -R openkim /home/openkim
USER openkim
ENV SHELL=/bin/bash
