### Publishing

We are currently migrating from using `publishDir` to workflow outputs.

- Workflows and subworkflows need to emit all result files that need to be published.
- Main workflow has logic to mix together files that will be published in the same output directory.
- Publish block should publish one channel per directory whenever possible in the main workflow (`main.nf`).
- Output block determines directories for each channel that comes from the publish block (`main.nf`).
- Keep file naming/prefixes in the configs as much as possible.
- Remove the corresponding `publishDir` entry from `conf/modules/` when adding a process to `ch_publish`.
