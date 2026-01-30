FROM continuumio/miniconda3:latest

COPY modules/local/gene_onthology/environment.yaml /tmp/environment.yaml

RUN conda env create -f /tmp/environment.yaml -n go_term_finder && \
    conda clean -afy

ENV PATH=/opt/conda/envs/go_term_finder/bin:$PATH

WORKDIR /work
CMD ["/bin/bash"]
