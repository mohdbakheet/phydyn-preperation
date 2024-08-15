# Using phydynBEAST for Multi-Host SIR Models

This README provides instructions on how to use the `phydynBEAST` package developed by Erik Volz to generate BEAST2 XML files for a Susceptible-Infected-Recovered (SIR) compartmental model across multiple host species.
We consider the following model

$$
\frac{dS_C}{dt} = - \left( \beta_CC \frac{S_C I_C}{N_C} + \beta_CC \frac{S_C I_C}{N_C}  + \beta_CC \frac{S_C I_C}{N_C} \right)
\frac{dE_C}{dt} = \left( \beta_CC \frac{S_C I_C}{N_C} + \beta_CC \frac{S_C I_C}{N_C}  + \beta_CC \frac{S_C I_C}{N_C} \right) - \gamma_C E_C

$$

## Installation

### 1. Install R
Ensure that R is installed on your computer. You can download it from [CRAN](https://cran.r-project.org).

### 2. Install phydynBEAST
Install the `phydynBEAST` package using `devtools`:
```R
if (!require("devtools")) install.packages("devtools")
devtools::install_github('emvolz/phydynbeast')
```

## Setting Up The Model

### 1. Load the Package
Load `phydynBEAST` into your R session:
```R
library(phydynBEAST)
```

### 2. Define the Multi-Host SIR Model
Create the compartmental model using epidemiological equations. Here is an example setup for a multi-host system:
```r
eqns <- list(
    confeqn('N_C = S_C + E_C + C + R_C', type='definition'),
  confeqn('N_T = S_T + E_T + T + R_T', type='definition'),
  confeqn('N_H = S_H + E_H + H + R_H', type='definition'),
 
  # Infection dynamics within and between species
  confeqn('beta_CC * C / N_C', origin='C', type='birth', destination='E_C'),
  confeqn('beta_TC * T / N_T', origin='T', type='birth', destination='E_T'),
  confeqn('beta_HC * H / N_H', origin='H', type='birth', destination='E_H'),
  
  
  confeqn('beta_CT * C / N_C', origin='C', type='birth', destination='E_db'),
  confeqn('beta_TT * T / N_T', origin='T', type='birth', destination='E_db'),
  confeqn('beta_HT * H / N_H', origin='H', type='birth', destination='E_db'),
 
  confeqn('beta_CH * C / N_C', origin='C', type='birth', destination='E_wm'),
  confeqn('beta_TH * T / N_T', origin='T', type='birth', destination='E_wm'),
  confeqn('beta_HH * H / N_H', origin='H', type='birth', destination='E_wm'),
  
  # Recovery and mortality
  
  confeqn('(gamma_C + mu_H) * C',  origin='C', type='death'),
  confeqn('(gamma_T + mu_H) * T',  origin='T', type='death'),
  confeqn('(gamma_H + mu_H) * H',  origin='H', type='death'),
  
  confeqn('eta_C * E_T', origin='E_C', type='migration', destination='C'),
  confeqn('eta_T * E_T', origin='E_T', type='migration', destination='T'),
  confeqn('eta_H * E_H', origin='E_H',  type='migration', destination='H'),
  
  confeqn('gamma_C * C', origin='R_C', type='nondeme'),
  confeqn('gamma_T * T', origin='R_T', type='nondeme'),
  confeqn('gamma_H * H', origin='R_H',  type='nondeme' ),
  
  # Susceptibility
  confeqn('-(beta_CC * C / N_C + beta_TC * T / N_T + beta_HC * H / N_H) * S_C', origin='S_C', type='nondeme'),
  confeqn('-(beta_CT * C / N_C + beta_TT * T / N_T + beta_HT * H / N_H) * S_T', origin='S_db', type='nondeme'),
  confeqn('-(beta_CH * C / N_C + beta_TH * T / N_T + beta_HH * H / N_H) * S_H', origin='S_h', type='nondeme')
)
```

### 3. Set Initial Conditions and Parameters
Define the initial conditions for each compartment and the parameter values:
```R
# Initial conditions
#init_conditions <- c(S_C = 1000, I_C = 10, R_C = 0)

# Parameters
#param_values <- c(beta_CC = 0.1, gamma_C = 0.05)
```

## Generating the BEAST2 XML File

### 1. Generate XML
Generate the XML configuration for the BEAST2 analysis:
```R
#xml <- generateXML(eqns, init_conditions, param_values, endTime = 50)
```

### 2. Save the XML to a File
Save the generated XML to a file for use with BEAST2:
```R
#write(xml, file = "SIR_MultiHost.xml")
```

## Running the Analysis

Once the XML file is prepared, you can run the analysis using BEAST2.

### 1. Run BEAST2
- Open BEAST2 and select the `SIR_MultiHost.xml` file.
- Configure your run settings according to your research needs.
- Start the analysis.
