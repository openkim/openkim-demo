# base image
FROM ghcr.io/openkim/developer-platform@sha256:aaef53d490b5f69798874ae0d929f2e5dfb1ddc8ba053d499256dbedce4d08c1

USER root
RUN pip install --no-cache-dir notebook jupyterlab

# Make sure the contents of our repo are in ${HOME}
COPY . /home/openkim/openkim-demo/
RUN chown -R openkim /home/openkim
USER openkim
ENV SHELL=/bin/bash
ENV LD_LIBRARY_PATH=:/usr/local/lib
