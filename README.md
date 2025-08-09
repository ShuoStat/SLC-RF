# Supervised Linear Combination Random Forest

We propose a novel random forest-based method with a two-layer structure, called the Supervised Linear Combination Random Forest (SLC-RF), which achieves impressive performance in high-dimensional data. In contrast to standard approaches, SLC-RF first transforms the input data into supervised linear combinations (Layer 1), improving feature representation based on the target variable. These enriched features are then fed into a conventional random forest (Layer 2) for prediction. The core innovation lies in simultaneously strengthening individual trees through supervised feature construction and promoting tree diversity via randomization applied at both layers.

### Project Structure

```
.
├── codes/
├── data/
├── output/
├── results/
└── README.md
```
`codes` Scripts for reproducing the results. Inside are codes for both real data application and simulation study.      
`data` Real-world datasets used for model development.    
`output` Intermediate results and temporary outputs.      
`results` Final processed results. 


### Reproducing the Results

Run *1_runReal.R* step by step to generate the results of real-world application  
Run *2_sumReal.R* to summarize the results into tables and figures
```
R /codes/realData/1_runReal.R
R /codes/realData/2_sumReal.R
```
Run *1_runSim.R* step by step to generate the results of simulations  
Run *2_sumSim.R* to summarize the simulation results into tables and figures 

```
R /codes/simulation/1_runSim.R
R /codes/simulation/2_sumSim.R
```

Of note, run on R version 4.3.3. Results may be slightly different with different R version.      

### Contact

Maintainer: Shuo Wang  
GitHub: https://github.com/ShuoStat/  
Email: wangsures@foxmail.com  

### License

### Custom License

BY-NC-ND 4.0 License © 2025 Shuo Wang

This work is provided for personal and educational use only.

You may not:
- Use this work in any form of academic publication, including but not limited to journal articles, conference papers, or theses.
- Use this work for any commercial purposes, including but not limited to products, services, or monetized platforms.
- Distribute modified or derivative versions of this work.

All rights are reserved. For inquiries, contact wangsures@foxmail.com

