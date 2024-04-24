Describe how to configure the workflow (using config.yaml and maybe additional files).
All of them need to be present with example entries inside of the config folder.


# General settings
To configure this workflow, modify ``config/config.yaml`` according to your needs, following the explanations provided in the file.

# Sample sheet

Add samples to `config/samples.tsv`. For each sample, the columns `sample_name`, `comparisons` have to be defined. 
* If multiple sample comparisons are needed, the same has to be assign to the particular columns
* No cell values has to be left blank, in case if comparison is between only two samples. In such cases, the blank cells has to be filled using `skip` word.
