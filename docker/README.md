This folder contains Docker build instructions for local modules used by the pipeline.

Build commands (run from repo root):

# Parser (uses modules/local/parser/envirionment.yaml)
```bash
docker build -f modules/local/parser/Dockerfile -t dogay/gb-parser:3.11 .
```

# GO term finder
```bash
docker build -f modules/local/gene_onthology/Dockerfile -t dogay/go-term-finder:3.11 .
```

# KEGG requests
```bash
docker build -f modules/local/kegg/Dockerfile -t dogay/kegg-requests:3.11 .
```

# Uniprot mapping
```bash
docker build -f modules/local/uniprot/Dockerfile -t dogay/uniprot-mapping:3.11 .
```

After building, push to your registry:
```bash
# example push (replace with your registry creds)
docker push dogay/gb-parser:3.11
# etc.
```

Then update the `container` fields in the corresponding `modules/local/*/main.nf` files to use the built images, for example:
```
container 'docker.io/dogay/gb-parser:3.11'
```

Notes:
- These Dockerfiles create conda environments inside the image from the module `environment.yaml` files. They are suitable for the `docker` profile.
- For smaller images consider using `mamba` or `micromamba` based base images.
