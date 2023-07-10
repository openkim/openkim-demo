# base image
FROM FROM ghcr.io/ilia-nikiforov-umn/developer-platform:1.2.2-ilia-dev-labels-0.1-minimal

USER root
RUN pip install --no-cache-dir notebook jupyterlab

# Make sure the contents of our repo are in ${HOME}
COPY . /home/openkim/openkim-demo/
RUN chown -R openkim /home/openkim
USER openkim
ENV SHELL=/bin/bash
ENV LD_LIBRARY_PATH=:/usr/local/lib
