# base image
FROM ghcr.io/ilia-nikiforov-umn/developer-platform@sha256:545d348368f79d9416f201a35b2c80ac365a205c4ad8614808a9c657fc78c55c

USER root
RUN pip install --no-cache-dir notebook jupyterlab

# Make sure the contents of our repo are in ${HOME}
COPY . /home/openkim/openkim-demo/
RUN chown -R openkim /home/openkim
USER openkim
ENV SHELL=/bin/bash
ENV LD_LIBRARY_PATH=:/usr/local/lib
