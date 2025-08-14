# base image
FROM ghcr.io/openkim/developer-platform@sha256:b1f2cee013e9d59aa1643cc8219115e8a9d01d1ee6fef99f28bd3bd018cd10d7

USER root
RUN pip install --no-cache-dir notebook jupyterlab==4.3.6 nglview==3.1.2

# Make sure the contents of our repo are in ${HOME}
COPY . /home/openkim/openkim-demo/
RUN chown -R openkim /home/openkim
USER openkim
ENV SHELL=/bin/bash
