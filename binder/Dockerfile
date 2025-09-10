# base image
FROM ghcr.io/openkim/developer-platform@sha256:ed19d309eda4f51d9f19d66fee2c37354e703e269d5dba3c3359f6bce3d6ffdb

USER root
RUN pip install --no-cache-dir notebook jupyterlab==4.3.6 nglview==3.1.2

# Make sure the contents of our repo are in ${HOME}
COPY . /home/openkim/openkim-demo/
RUN chown -R openkim /home/openkim
USER openkim
ENV SHELL=/bin/bash
