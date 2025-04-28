# base image
FROM ghcr.io/openkim/developer-platform@sha256:925802c5cac32b1e9b0f7770c633be8094a11fe7af467c4fd49b7fc4766c5954

USER root
RUN pip install --no-cache-dir notebook jupyterlab==4.3.6 nglview==3.1.2

# Make sure the contents of our repo are in ${HOME}
COPY . /home/openkim/openkim-demo/
RUN chown -R openkim /home/openkim
USER openkim
ENV SHELL=/bin/bash
