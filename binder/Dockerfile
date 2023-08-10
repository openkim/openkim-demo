# base image
FROM ghcr.io/ilia-nikiforov-umn/developer-platform@sha256:b8a591ca1d027b7761353a062b6a21c27104fc65c49fe5eb50ca0fa308cbc55b

USER root
RUN pip install --no-cache-dir notebook jupyterlab nglview

# Make sure the contents of our repo are in ${HOME}
COPY . /home/openkim/openkim-demo/
RUN chown -R openkim /home/openkim
USER openkim
ENV SHELL=/bin/bash
ENV LD_LIBRARY_PATH=:/usr/local/lib
