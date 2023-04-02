# base image
FROM ghcr.io/openkim/developer-platform:latest-minimal

USER root
RUN apt-get update && apt-get install -y git
RUN pip install --no-cache-dir notebook jupyterlab
RUN pip install requests
RUN cd /home/openkim && git clone https://github.com/ilia-nikiforov-umn/mach-2023-openkim-tutorial
