# `genomic-medicine-sweden/nallo`: Contributing Guidelines

Hi there!
Many thanks for taking an interest in improving genomic-medicine-sweden/nallo.

We try to manage the required tasks for genomic-medicine-sweden/nallo using GitHub issues, you probably came to this page when creating one.
Please use the pre-filled template to save time.

However, don't be put off by this template - other more general issues and suggestions are welcome!
Contributions to the code are even more welcome ;)

## Table of contents

- [General](#general)
  - [Contribution workflow](#contribution-workflow)
    - [Pull Requests](#pull-requests)
      - [Review](#review)
  - [Software versioning, changelog and updates](#software-versioning-changelog-and-updates)
    - [Semantic versioning and changelog](#semantic-versioning-and-changelog)
      - [Patch](#patch)
    - [Nextflow version bumping](#nextflow-version-bumping)
    - [Update nf-core template](#update-nf-core-template)
  - [Developer setup](#developer-setup)
    - [Installation and dependencies for development](#installation-and-dependencies-for-development)
    - [GitHub Codespaces](#github-codespaces)
  - [Running tests](#running-tests)
    - [Lint tests](#lint-tests)
    - [Pipeline tests](#pipeline-tests)
  - [Adding new tools](#adding-new-tools)
    - [1. Update citations](#1-update-citations)
    - [2. Update `README.md`](#2-update-readmemd)
- [Coding conventions](#coding-conventions)
  - [Architecture & structure](#architecture--structure)
  - [Adding a new step](#adding-a-new-step)
  - [Channels](#channels)
    - [Naming schemes](#naming-schemes)
  - [Parameters](#parameters)
  - [Publishing](#publishing)
  - [Configuration](#configuration)
  - [Writing tests](#writing-tests)
  - [Style](#style)

## General

### Contribution workflow

If you'd like to write some code for genomic-medicine-sweden/nallo, the standard workflow is as follows:

1. Check that there isn't already an issue about your idea in the [genomic-medicine-sweden/nallo issues](https://github.com/genomic-medicine-sweden/nallo/issues) to avoid duplicating work. If there isn't one already, please create one so that others know you're working on this
2. [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [genomic-medicine-sweden/nallo repository](https://github.com/genomic-medicine-sweden/nallo) to your GitHub account
3. Make the necessary changes / additions within your forked repository following [Pipeline conventions](#pipeline-contribution-conventions)
4. Use `nf-core pipelines schema build` and add any new parameters to the pipeline JSON schema (requires [nf-core tools](https://github.com/nf-core/tools) >= 1.10).
5. Submit a Pull Request against the `dev` branch and wait for the code to be reviewed and merged

If you're not used to this workflow with git, you can start with some [docs from GitHub](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests) or even their [excellent `git` resources](https://try.github.io/).

#### Pull Requests

When opening a pull request to suggest changes to the code, please make sure to follow the [Pipeline contribution conventions](#pipeline-contribution-conventions) for the code and to fill in the necessary information in the pull request template as well as address all points in the `PR checklist`.

##### Review

When reviewing a PR, make sure to check that:

- The code follows the [Pipeline contribution conventions](#pipeline-contribution-conventions).
- The information in the PR (and related issue) is clear and sufficient to understand the change and the motivation for it - title, description and entry in `CHANGELOG.md`, if applicable.
- All the items in the `PR checklist` have been addressed, the changes are well documented and the tests are passing.

Be positive and constructive in your review, and whenever possible offer suggestions for improvement rather than just pointing out issues.

### Software versioning, changelog and updates

#### Semantic versioning and changelog

Release versioning tries to be maintained according to [Semantic Versioning](https://semver.org/spec/v2.0.0.html) and a changelog is maintained according to the [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) format. Altough, for a Nextflow pipeline it's hard do decide what is a breaking/non-breaking API change.

The `Fixed` section of the changelog should be reserved for bugs fixed from one release to another, but bug fixes that appeared and was solved in dev. These should go in the `Changed` section.

##### Patch

:warning: Only in the unlikely and regretful event of a release happening with a bug.

- On the genomic-medicine-sweden/nallo repository, make a new branch `patch` based on `upstream/main` or `upstream/master`.
- Fix the bug, and bump version (X.Y.Z+1).
- Open a pull-request from `patch` to `main`/`master` with the changes.

#### Nextflow version bumping

You may bump the minimum required version of nextflow in the pipeline with: `nf-core pipelines bump-version --nextflow . [min-nf-version]`

#### Update nf-core template

Since this is not an nf-core pipeline, the nf-core template is not automatically updated in the `TEMPLATE` branch. Follow these step to update the template:

1. Update the `TEMPLATE` branch by running `nf-core pipelines sync`. Fix any merge conflicts and open a PR to then merge the changes.
1. Open a PR to merge the `TEMPLATE` branch into `dev` to update the template files in the main codebase.

### Developer setup

#### Installation and dependencies for development

In order to run the pipeline, develop and test your changes locally, we recommend that you set up:

- A conda environment with `nextflow`, `nf-core` tools and `nf-test`. - For this, follow the instructions from the [nf-core documentation](https://nf-co.re/docs/nf-core-tools/cli/installation#install-with-conda) and install `nf-test` from bioconda or by following the [nf-test installation instructions](https://www.nf-test.com/installation/) . Additional information about [installation of nf-core dependencies](https://nf-co.re/docs/usage/getting_started/installation/) is also available, if needed.
- Install Docker (https://www.docker.com/products/docker-desktop/) and make sure the daemon is running when you want to run the tests locally.

Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

Additionally, pre-commit hooks are set up to automatically check the code and generate parameters documentation when committing. To install the pre-commit hooks, run `pre-commit install` in the root of the repository. Note that, other than the default pre-commit hooks from the nf-core template, there is an additional hook - [`nf-core pipelines schema docs` pre-commit hook](https://github.com/genomic-medicine-sweden/nf-core-schema-docs) - set up in this codebase for automatically generating parameters documentation if there are any changes to parameters in your commit.

#### GitHub Codespaces

This repo includes a devcontainer configuration which will create a GitHub Codespaces for Nextflow development! This is an online developer environment that runs in your browser, complete with VSCode and a terminal.

To get started:

- Open the repo in [Codespaces](https://github.com/genomic-medicine-sweden/nallo/codespaces)
- Tools installed
  - nf-core
  - Nextflow

Devcontainer specs:

- [DevContainer config](.devcontainer/devcontainer.json)

### Running tests

You have the option to test your changes locally by running the pipeline. For receiving warnings about process selectors and other `debug` information, it is recommended to use the debug profile. Execute all the tests with the following command:

```bash
nf-test test --profile debug,test,docker --verbose
```

You can also run `nf-test test <path>` for a single given test, `nf-test test --tag PIPELINE` to run all pipeline tests, or `nf-test test --tag call_snvs,annotate_snvs` to run multiple subworkflow test.

When you create a pull request with changes, [GitHub Actions](https://github.com/features/actions) will run automatic tests.
Typically, pull-requests are only fully reviewed when these tests are passing, though of course we can help out before then.

There are typically two types of tests that run:

#### Lint tests

`nf-core` has a [set of guidelines](https://nf-co.re/developers/guidelines) which all pipelines must adhere to.
To enforce these and ensure that all pipelines stay in sync, they have developed a helper tool which runs checks on the pipeline code. This is in the [nf-core/tools repository](https://github.com/nf-core/tools) and once installed can be run locally with the `nf-core pipelines lint <pipeline-directory>` command.

If any failures or warnings are encountered, please follow the listed URL for more documentation.

#### Pipeline tests

This pipeline is set up with a minimal set of test-data.
`GitHub Actions` then runs the pipeline on this data to ensure that it exits successfully.
If there are any test failures then the automated check has status set to fail.
These tests are run both with the latest available version of `Nextflow` and also the minimum required version that is stated in the pipeline code.

### Adding new tools

#### 1. Update citations

When adding a new tool to the pipeline, update `assets/software_references.yml` with the citation and bibliography.

#### 2. Update `README.md`

Add the tool to the relevant numbered section in the **Pipeline summary**. If the tool belongs to a new category not yet represented, add a new numbered section in the appropriate position.

## Coding conventions

To make the `genomic-medicine-sweden/nallo` code and processing logic more understandable for new contributors and to ensure quality, we semi-standardise the way the code and other contributions are written.

### Architecture & structure

- **Use subworkflows** — Don't add logic to `workflows/nallo.nf` that is specific to a subworkflow. Create new subworkflows as needed under `subworkflows/` and import them into `workflows/nallo.nf`.
- **Reuse over duplication** — Follow the DRY (Don't Repeat Yourself) principle.
- **nf-core modules take precedence** — Always prefer a module from nf-core or genomic-medicine-sweden over writing a local one. Only add to modules/local/ as a last resort when the use case is too pipeline-specific.
- **Use and share subworkflows with the GMS community** — subworkflows from [genomic-medicine-sweden/nf-core-modules](https://github.com/genomic-medicine-sweden/nf-core-modules) are intended for use across pipelines within the Genomic Medicine Sweden group. Prefer using and contributing to these rather than writing pipeline-specific code in `subworkflows/local/`. If you think a subworkflow could be useful for other pipelines, consider adding it there instead of `subworkflows/local/`.

### Adding a new step

If you wish to contribute a new step, please use the following coding standards:

1. Define the corresponding input channel into your new process from the expected previous process channel.
2. Write the process block (see below).
3. Define the output channel if needed (see below).
4. Add any new parameters to `nextflow.config` with a default (see below).
5. Add any new parameters to `nextflow_schema.json` with help text (via the `nf-core pipelines schema build` tool).
6. Add sanity checks and validation for all relevant parameters.
7. Perform local tests to validate that the new code works as expected.
8. If applicable, add a new test in the `tests` directory.
9. Update MultiQC config `assets/multiqc_config.yml` so relevant suffixes, file name clean up and module plots are in the appropriate order. If applicable, add a [MultiQC](https://https://multiqc.info/) module.
10. Add a description of the output files and if relevant any appropriate images from the MultiQC report to `docs/output.md`.

### Channels

- **Conditional channels**: always initialize to `channel.empty()` before any `if` block that may or may not assign them. Never leave a channel potentially undefined.

### Parameters

- `params` must only be accessed in the main unnamed workflow (`workflow` in `main.nf`). Subworkflows and named workflows receive all values as explicit `val_*` arguments. Never reference `params` directly inside a subworkflow.
- Parameters should be initialised/defined with default values within the `params` scope in `nextflow.config`. Don't hardcode values that a user might reasonably want to change. Once added, run nf-core pipelines schema build to register them in nextflow_schema.json.

<!-- TODO: Add information about parameters for skipping tools when that logic is decided. -->

### Publishing

We are currently migrating from using `publishDir` to workflow outputs.

- Workflows and subworkflows needs to emit all result files that need to be published.
- Main workflow has logic to mix together files that will be published in the same output directory.
- Publish block should publish one channel per directory whenever possible in the main workflow (main.nf)
- Output block determine directories for each channel that comes from publish block (main.nf).
- Keep file naming/prefixes in the configs as much as possible.
- Remove the corresponding `publishDir` entry from `conf/modules/` when adding a process to `ch_publish`.

### Configuration

- Process-level options go in `conf/subworkflows/<subworkflow_name>.config`, not inline in the subworkflow `.nf` file.
- Use module configs strictly for defining `ext.args`, `ext.args2`, and `ext.prefix`.
- Conditional behavior (e.g. save as CRAM vs BAM) should be handled in the subworkflow — not via config-level flags.
- Process resource requirements (CPUs / memory / time) go in `conf/base.config` using `withLabel:` selectors so they can be shared across processes. Use `${task.cpus}` and `${task.memory}` in `script:` blocks to apply them dynamically.

### Writing tests

- Every subworkflow should have real tests and `-stub` tests at `subworkflows/local/<name>/tests/main.nf.test`.
- Snapshot files (`*.nf.test.snap`) are committed alongside tests — update them when outputs change.
- Pipeline-level tests live in `tests/`.
- When modules in tests require different parameters, use `params { <module>_args = ... }` similar to https://nf-co.re/docs/specifications/components/modules/testing#configuration-of-extargs-in-tests instead of having mutliple configs.

### Style

- Both `take:` and `emit:` block entries require an inline type comment. Use `name // type: [mandatory|optional] description` for `take:` and `name = value // channel: [type description]` for `emit:`. Always include the comment — never leave an entry uncommented.

  ```groovy
  take:
  ch_vcf                // channel: [mandatory] [ val(meta), path(vcf) ]
  ch_reduced_penetrance // channel: [optional]  [ path(penetrance) ]
  val_aligner           // string:  [mandatory] aligner name (bwa/bwamem2/bwameme)
  process_with_sort     // Boolean

  emit:
  vcf     = ch_vcf      // channel: [ val(meta), path(vcf) ]
  publish = ch_publish  // channel: [ val(destination), val(value) ]
  ```

- Avoid using `ch_* = <...>`, use the `.set {ch_*}` operator to create new channels whenever possible.

- Use the nextflow code formatting (`nextflow lint -harshil-alignment -format file.nf`) after making changes, but beware that it currently removes some inline comments.
- Use comment-blocks for multi-line comments, e.g.:

```
/*
 * This is
 * a multi-line comment.
 */
```

- Use `[...]` instead of `tuple(...)`.
- Avoid nested if-else statements, prefer single-level statements when possible.
