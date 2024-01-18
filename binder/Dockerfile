# base image
FROM ghcr.io/openkim/developer-platform@sha256:eaf31d5fbc01d16631c85683b2ee26d59990a1950113cb175bb3754c0f2f10fc

USER root
RUN pip install --no-cache-dir notebook jupyterlab nglview

# Make sure the contents of our repo are in ${HOME}
COPY . /home/openkim/openkim-demo/
RUN chown -R openkim /home/openkim
USER openkim
ENV SHELL=/bin/bash
ENV LD_LIBRARY_PATH=:/usr/local/lib
