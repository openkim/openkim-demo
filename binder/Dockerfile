# base image
FROM ghcr.io/openkim/developer-platform@sha256:ac55a97adbf8b5e844ba7f1a05560f07925e81383c828eb9148dd2937c73e637

USER root
RUN pip install --no-cache-dir notebook jupyterlab nglview

# Make sure the contents of our repo are in ${HOME}
COPY . /home/openkim/openkim-demo/
RUN chown -R openkim /home/openkim
USER openkim
ENV SHELL=/bin/bash
ENV LD_LIBRARY_PATH=:/usr/local/lib
