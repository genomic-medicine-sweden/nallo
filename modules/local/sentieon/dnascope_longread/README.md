## Container

This container was built with the following instructions:

```sh
docker build --build-arg SENTIEON_VERSION=202503.02 --build-arg SENTIEON_CLI_BRANCH=8bd23b6 -t clinicalgenomicslund/dnascope-longread:1.4.0-lrRPA-f9c8811 .
```

Where `SENTIEON_VERSION` is the version of the sentieon binary, and `SENTIEON_CLI_BRANCH` is the branch/release commit of the sentieon-cli.

The sentieon-cli package that contains the DNAscope longread snv calling pipeline used by the module. See here for more info: https://github.com/Sentieon/sentieon-cli/

## Dockerfile

The Dockerfile is based the Dockerfile supplied by sentieon here: https://github.com/Sentieon/sentieon-docker ([commit efb2a31](https://github.com/Sentieon/sentieon-docker/commit/efb2a31f3cdae86e15cc1007013cba3083858a1d))

It is extended with all dependencies that are required by sentieon-cli / dnascope-longread.
