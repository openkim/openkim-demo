# base image
FROM ghcr.io/openkim/developer-platform:latest-minimal

USER root
RUN apt-get update && apt-get install -y git
RUN pip install --no-cache-dir notebook jupyterlab
