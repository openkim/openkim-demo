# base image
FROM ghcr.io/openkim/developer-platform:latest-minimal

USER root
RUN pip install --no-cache-dir notebook jupyterlab

# Make sure the contents of our repo are in ${HOME}
COPY . /home/openkim/mach-2023-openkim-tutorial/
RUN chown -R openkim /home/openkim
USER openkim
ENV SHELL=/bin/bash
