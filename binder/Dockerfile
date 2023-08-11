# base image
FROM ghcr.io/openkim/developer-platform@sha256:bb1a32c0aa7147643f234be7b9341e4e984eea0de93f438ccc73459da7bc0e32

USER root
RUN pip install --no-cache-dir notebook jupyterlab nglview

# Make sure the contents of our repo are in ${HOME}
COPY . /home/openkim/openkim-demo/
RUN chown -R openkim /home/openkim
USER openkim
ENV SHELL=/bin/bash
ENV LD_LIBRARY_PATH=:/usr/local/lib
