FROM continuumio/miniconda3:latest

COPY modules/local/uniprot/environment.yaml /tmp/environment.yaml

RUN conda env create -f /tmp/environment.yaml -n uniprot_mapping && \
    conda clean -afy

ENV PATH=/opt/conda/envs/uniprot_mapping/bin:$PATH

WORKDIR /work
CMD ["/bin/bash"]
