# base image
FROM ghcr.io/ilia-nikiforov-umn/developer-platform@sha256:d7ca16a666b66f98e9f6a9f8a32ae7cb307674a01535737b495d8f1838e38628

USER root
RUN pip install --no-cache-dir notebook jupyterlab nglview

# Make sure the contents of our repo are in ${HOME}
COPY . /home/openkim/openkim-demo/
RUN chown -R openkim /home/openkim
USER openkim
ENV SHELL=/bin/bash
ENV LD_LIBRARY_PATH=:/usr/local/lib
