# base image
FROM ghcr.io/openkim/developer-platform@sha256:ab2cd848f8971930b88d89aa5459195e718a908895544f16e3613c7bd5efa4b5

USER root
RUN pip install --no-cache-dir notebook jupyterlab nglview

# Make sure the contents of our repo are in ${HOME}
COPY . /home/openkim/openkim-demo/
RUN chown -R openkim /home/openkim
USER openkim
ENV SHELL=/bin/bash
ENV LD_LIBRARY_PATH=:/usr/local/lib
