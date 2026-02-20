# TRENDYv12 SOC Change Analysis (Fork)
This repository is a fork of the TRENDYv12 S2 soil organic carbon (SOC) change analysis:
https://github.com/yinonbaron/SOC_change_TRENDYv12/tree/main/data

It was created to supplement the code required to reproduce results presented in:

> "A large global soil carbon sink informed by repeated soil samplings" 
> (unrevised preprint: https://www.biorxiv.org/content/10.1101/2025.04.25.650716v1).

# Modifications from the main branch
1. Simulation outputs are spatially masked to match the study domain used in the data-driven analysis described above.
2. Carbon pools of litter (`cLitter`) and woody debris (`cCwd`) are excluded. Only the soil pool (`cSoil`) is retained for more consistent comparison with the data-driven results.


