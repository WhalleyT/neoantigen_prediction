## Aims

## Installation
Make sure you have Nextflow installed. It can be downloaded [here](https://www.nextflow.io/).

Clone this repository by running `git clone https://github.com/whalleyt/neoantigen_prediction`. 
Assuming you are going to run this in one of the containers supplied, then you are ready to go.

## Usage
The parameters for this pipeline can be supplied at the command line or in nextflow.config. At it's most
basic one can run:

`nextflow main.nf`

There are a number of pre-supplied profiles to be used for various executor methods and containers. These can
be accesed by the `-profile` argument. For example to run on Slurm using a singularity container run:

`nextflow main.nf -profile slurm,singularity`

## Workflow

![workflow](assets/pipeline.png)

## References

## License

## Contact
