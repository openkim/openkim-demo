# base image
FROM ghcr.io/openkim/developer-platform@sha256:5ec388d5b7201a2bf8c29d1f2730a6145ba7e540b30b21e64c3365d0b964606b

USER root
RUN pip install --no-cache-dir notebook jupyterlab==4.3.6 nglview==3.1.2

# Make sure the contents of our repo are in ${HOME}
COPY . /home/openkim/openkim-demo/
RUN chown -R openkim /home/openkim
USER openkim
ENV SHELL=/bin/bash
