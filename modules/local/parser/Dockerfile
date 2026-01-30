FROM continuumio/miniconda3:latest

# Copy the env file from the repo into the image
COPY modules/local/parser/envirionment.yaml /tmp/environment.yaml

# Create the conda env and clean up
RUN conda env create -f /tmp/environment.yaml -n gb_parser && \
    conda clean -afy

# Ensure the environment is on PATH
ENV PATH=/opt/conda/envs/gb_parser/bin:$PATH

WORKDIR /work
CMD ["/bin/bash"]
