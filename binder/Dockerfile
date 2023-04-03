# base image
FROM ghcr.io/openkim/developer-platform:latest-minimal

RUN pip install --no-cache-dir notebook jupyterlab

# Make sure the contents of our repo are in ${HOME}
COPY .. ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}


