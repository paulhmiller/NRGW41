NRG-W41 Project
===============

Project lead: Paul Miller [paulhmiller@gmail.com, pmiller@bccrc.ca]

## Summary

C57BL/6-W41 mice were crossed with NRG mice, the het progeny intercrossed, and
then mice homozygous for NOD^Sirp, Rag1-/-, IL2Rg-/-, and Kit^W41 were
selected. For the tranplants done by Paul, no further backcrossing was
performed (although 2 backcrossing projects, to NRG and C57Bl/6, were done for
future experiments). The former (mixed background / "dirty") mice were also
crossed with NRG-3GS mice to generate NRG-W41-3GS mice, which Naoto has used
as recipients for leukemic cells. But this current project is for the NRG-W41 mice,
transplanted with CB cells. There are also BM transplanted mice but there
is only one replicate.

Variables include:
- Strain (NRG versus NRG-W41)
- Irradiated (no irradiation versus 900 cGy NRG versus 150 cGy NRG-W41)
- Sex (both strains)

The cell dose used was either 2 or 4 x 10e4 CD34 per mouse. Some experiments
used barcoded cells. Note that the 40K cell dose was actually 38K and is
sometimes written so. 


## Notation

Experiment names contain a date in YYMMDD format with a letter suffix.
Sometimes only the letter is used. But as I cycled through A-Z for my PhD
there are muliple experiments for each letter. Difficult to mix them up, but
still important to know. These experiment codes are important to distinguish
replicates and also to go back to the raw data. 


## Sub-directory Information

./figures
- contains plotting R scripts and plots. T

./primary_kinetics_pool
- data for primary transplants. CSV files each contain data from multiple 
  experiments, and have both 40K CB CD34+ data and (pooled) 20K and 40K CB data. 
  My analysis shows that the results for 20K and 40K are similar, so pooling 
  is appropriate. There is also data from a single BM transplant experiment. 
  All of these have both BM and PB data sets. 

./LDA_secondary
- data and LDA calculation scripts from secondary transplants.

./barcoding
- information on the barcoded samples. Currently a bit of a mess as I am still
  working on this. 

./reports
- text and presentation files for external reporting (meetings, grants, etc). 


## Miscellaneous Information
- Units of PB data are in microlitres
- Paul's preferred PB detection threshold is 0.2/uL for leukocytes and 25/uL
  for platlets. 

