# base image
FROM ghcr.io/openkim/developer-platform:latest-minimal

USER root
RUN pip install --no-cache-dir notebook jupyterlab

# Make sure the contents of our repo are in ${HOME}
COPY . ${HOME}
RUN chown openkim ${HOME}
USER openkim


