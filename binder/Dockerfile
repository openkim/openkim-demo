# base image
FROM ghcr.io/openkim/developer-platform@sha256:b5f63972d307e74ccc03e237a6eb3ca51502a44fb34a96999fe6c562cfe73e4c

USER root
RUN pip install --no-cache-dir notebook jupyterlab nglview==3.1.2

# Make sure the contents of our repo are in ${HOME}
COPY . /home/openkim/openkim-demo/
RUN chown -R openkim /home/openkim
USER openkim
ENV SHELL=/bin/bash
