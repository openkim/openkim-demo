# base image
FROM ghcr.io/openkim/developer-platform@sha256:20e579c9752fee272c60778183dd0ab0d0de1dd04a8a09952bd59993d8b12c78

USER root
RUN pip install --no-cache-dir notebook jupyterlab==4.3.6 nglview==3.1.2

# Make sure the contents of our repo are in ${HOME}
COPY . /home/openkim/openkim-demo/
RUN chown -R openkim /home/openkim
USER openkim
ENV SHELL=/bin/bash
