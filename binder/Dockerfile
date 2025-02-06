# base image
FROM ghcr.io/openkim/developer-platform@sha256:04b9c0baa692effc55a527a3fad0eb671c74ee1572ffca89a1ca642c6b5c7f5d

USER root
RUN pip install --no-cache-dir notebook jupyterlab nglview==3.1.2

# Make sure the contents of our repo are in ${HOME}
COPY . /home/openkim/openkim-demo/
RUN chown -R openkim /home/openkim
USER openkim
ENV SHELL=/bin/bash
