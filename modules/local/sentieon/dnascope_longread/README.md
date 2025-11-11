## Container

This container was built with the following instructions:

```sh
docker build --build-arg SENTIEON_VERSION=202503.01 --build-arg SENTIEON_CLI_BRANCH=f9c8811 -t clinicalgenomicslund/dnascope-longread:1.4.0-lrRPA-f9c8811 .
```

Where `SENTIEON_VERSION` is the version of the sentieon binary, and `SENTIEON_CLI_BRANCH` is the branch/commit of the sentieon-cli package that includes the DNAscope Longread snv caller. See here for more info: https://github.com/Sentieon/sentieon-cli/

## Dockerfile

The Dockerfile is based the Dockerfile supplied by sentieon here: https://github.com/Sentieon/sentieon-docker ([commit efb2a31](https://github.com/Sentieon/sentieon-docker/commit/efb2a31f3cdae86e15cc1007013cba3083858a1d))

It is extended with all dependencies that are required by sentieon-cli / dnascope-longread.
