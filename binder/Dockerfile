# base image
FROM ghcr.io/openkim/developer-platform:latest-minimal

USER root
RUN pip install --no-cache-dir notebook jupyterlab

# Make sure the contents of our repo are in ${HOME}
ADD . ${HOME}/mach-2023-openkim-tutorial
RUN chown -R openkim ${HOME}
USER openkim
