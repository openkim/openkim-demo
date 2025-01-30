# base image
FROM ghcr.io/openkim/developer-platform@sha256:720df3ce8024c46edaf1fda37d6ac26822a29cf9952bad8080041688721eb0cb

USER root
RUN pip install --no-cache-dir notebook jupyterlab nglview

# Make sure the contents of our repo are in ${HOME}
COPY . /home/openkim/openkim-demo/
RUN chown -R openkim /home/openkim
USER openkim
ENV SHELL=/bin/bash
ENV LD_LIBRARY_PATH=:/usr/local/lib
