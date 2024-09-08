# Assessing early detection ability through spatial arrangements in environmental surveillance for poliovirus: a simulation-based study

This repository contains the data and code for our paper:
> Asakura, T. R., & O'Reilly, K. M. (2024). Assessing early detection ability through spatial arrangements in environmental surveillance for poliovirus: a simulation-based study. medRxiv, 2024-08.

### Data used in our study
- Proportion of individuals shedding any amount of virus, [Tebbens RJD 2013](https://doi.org/10.1111/risa.12031):
`/data/Radboud_2013_PT.csv`
- Wastewater testing data, [Fuqing Wu 2021](https://doi.org/10.1016/j.watres.2021.117400):
`/data/fuqing_wu_2021.csv`

- OPV vaccine coverage in South Africa
: `/data/zaf_OPV_HEXA_vaccine_coverage_2020.csv`

- Wastewater plant information
: `/data/ES_surveillance_information_20240318.xlsx`

- Merged population data from [WolrdPop](https://www.worldpop.org/)
:`/data_pop/*.tif`
    - Code for merging origianl data is located in `2_1j` file.

- Boundary data from [geoBoundary](https://www.geoboundaries.org/):
`/dt_geoBoundaries-*`

- Simulation results to reproduce the figures: `/dt_tmp_hpc`

### How to run the code
Clone this repository and type `docker compose up` to
install the Docker image and set up the Docker container.
Then, you can run the code in `src` via Jupyter Lab.

Run the code in a sequence following the prefix number of each file in `src` directory.
It takes around 5-8 hours for each scenario to be completed for `4j` files
(which means several days are required to complete all the simulations).
It is recommended that the parallel computing method is implemented. See `slurm_v2` branch for example.

### License
[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE)