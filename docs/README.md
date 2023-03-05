# Documentation

## Example `samplesheet.csv`

```
group,input_file
reference,SAMPLE01.bam
reference,SAMPLE02.bam
query,SAMPLE03.bam
```

## Installation

### Bianca (UPPMAX)

1. For usage on Bianca (UPPMAX) make sure you have [`nextflow`](https://www.nextflow.io) in your path.
2. Have your wharf mounted on a computer where you have internet access (see the [transit-user-guide](https://www.uppmax.uu.se/support/user-guides/transit-user-guide/))
3. Change the nextflow storage location to somewhere on your mounted wharf, i.e. by running ```export NXF_ASSETS=/path/to/wharf/assets/```
4. Fetch the pipeline by running: ```nextflow pull fellen31/skierfe -r dev```
5. Pull required containers and set `NXF_SINGULARITY_CACHEDIR` or `NXF_SINGULARITY_LIBRARYDIR`? 
6. Log in to BIANCA and move the contents of `assets` from your wharf to where you want to store the pipeline. `NXF_ASSETS` points to `$NXF_HOME/assets` by default, make sure to set and export the `NXF_ASSETS` variable to the correct path if you choose to store them somewhere else.
7. You can now move to your project directory and run the pipeline by executing ```nextflow run fellen31/skierfe -r dev``` 

## Usage

Run the pipeline with: ```nextflow run fellen31/skierfe -r dev``` 

