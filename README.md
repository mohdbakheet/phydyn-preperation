# Using phydynBEAST for Multi-Host SIR Models

This README provides instructions on how to use the `phydynBEAST` package developed by Erik Volz to generate BEAST2 XML files for a Susceptible-Infected-Recovered (SIR) compartmental model across multiple host species.
We consider the following model

$$
\frac{dS_C}{dt} = - \left( \beta_CC \frac{S_C I_C}{N_C} + \beta_CC \frac{S_C I_C}{N_C}  + \beta_CC \frac{S_C I_C}{N_C} \right) \frac{-b \pm \sqrt{b^2 - 4ac}}{2a}
$$


$$
 = 
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

```r
parms <- list(
  confparm('beta_TC',
           initial = 3,
           prior = 'lognormal',
           operator = 'realrw',
           lower = 1,
           upper = 10,
           M = 1.09861229, 
           S = 0.5
  ),
  confparm('beta_CC',
           initial = 2.5,
           prior = 'lognormal',
           operator = 'realrw',
           lower = 1,
           upper = 10,
           M = 0.91629073, 
           S = 0.5
  ),
  confparm('beta_HC',
           initial = 4,
           prior = 'lognormal',
           operator = 'realrw',
           lower = 1,
           upper = 10,
           M = 1.38629436,
           S = 0.5
  ),
  confparm('beta_CT',
           initial = 4.5,
           prior = 'lognormal',
           operator = 'realrw',
           lower = 1,
           upper = 10,
           M = 1.5040774,
           S = 0.5
  ),
  confparm('beta_TT',
           initial = 5,
           prior = 'lognormal',
           operator = 'realrw',
           lower = 1,
           upper = 10,
           M = 1.60943791,
           S = 0.5
  ),
  confparm('beta_HT',
           initial = 3,
           prior = 'lognormal',
           operator = 'realrw',
           lower = 1,
           upper = 10,
           M = 1.09861229, 
           S = 0.5
  ),
  confparm('beta_CH',
           initial = 2.5,
           prior = 'lognormal',
           operator = 'realrw',
           lower = 1,
           upper = 10,
           M = 0.91629073, 
           S = 0.5
 ),
  confparm('beta_TH',
           initial = 4,
           prior = 'lognormal',
           operator = 'realrw',
           lower = 1,
           upper = 10,
           M = 1.38629436,
           S = 0.5
  ),
  confparm('beta_HH',
           initial = 4.5,
           prior = 'lognormal',
           operator = 'realrw',
           lower = 1,
           upper = 10,
           M = 1.5040774,
           S = 0.5
  ),
  confparm('eta_C',
           initial = 2,
           prior = 'lognormal',
           operator = 'realrw',
           lower = 1,
           upper = 3,
           M = 0.69314718,
           S = 0.2
  ),
  confparm('eta_T',
           initial = 2,
           prior = 'lognormal',
           operator = 'realrw',
           lower = 1,
           upper = 3,
           M = 0.69314718,
           S = 0.2
  ),
  confparm('eta_H',
           initial = 3,
           prior = 'lognormal',
           operator = 'realrw',
           lower = 2,
           upper = 4,
           M = 1.09861229,
           S = 0.2
  ),
 confparm('mu_C',
          initial = 0.75,
          prior = 'exponential',
          operator = 'realrw',
          lower = 0.5,
          upper = 1.0,
          mean = 0.75
 ),
 confparm('mu_T',
          initial = 0.95,
          prior = 'exponential',
          operator = 'realrw',
          lower = 0.9,
          upper = 1.0,
          mean = 0.95
 ),
 confparm('mu_H',
          initial = 0.75,
          prior = 'exponential',
          operator = 'realrw',
          lower = 0.5,
          upper = 1.0,
          mean = 0.75
 ),
  confparm('gamma_C',
           initial = (1 / 7),
           prior = 'lognormal',
           operator = 'realrw',
           lower =  (1 / 14),
           upper =  (1 / 3),
           M = -1.9459101490553135,
           S = 0.2
  ),
  confparm('gamma_T',
           initial = (1 / 5),
           prior = 'lognormal',
           operator = 'realrw',
           lower = (1 / 10),
           upper = (1 / 3),
           M = -1.6094379124341003,
           S = 0.2
  ),
  confparm('gamma_H',
           initial = (1 / 6),
           prior = 'lognormal',
           operator = 'realrw',
           lower = (1 / 7),
           upper = (1 / 5),
           M = log(365 * (1 / 6)),
           S = 0.2
  ),
  confparm('S_C',
           initial = 2e3,
           initial_condition_parameter = TRUE,
           prior = 'exponential',
           operator = 'realrw',
           estimate = TRUE,
           lower = 0,
           upper = 1e4,
           mean = 10
  ),
  confparm('E_C',
           initial = 1,
           initial_condition_parameter = TRUE,
           prior = 'exponential',
           operator = 'realrw',
           estimate = TRUE,
           lower = 0,
           upper = 1e4,
           mean = 10
  ),
  confparm('C',
           initial = 1,
           initial_condition_parameter = TRUE,
           prior = 'exponential',
           operator = 'realrw',
           estimate = TRUE,
           lower = 0,
           upper = 1e4,
           mean = 10
  ),
  confparm('R_C',
           initial = 0,
           initial_condition_parameter = TRUE,
           prior = 'exponential',
           operator = 'realrw',
           estimate = TRUE,
           lower = 0,
           upper = 1e3,
           mean = 10
  ),
  confparm('S_T',
           initial = 1500,
           initial_condition_parameter = TRUE,
           prior = 'exponential',
           operator = 'realrw',
           estimate = TRUE,
           lower = 0,
           upper = 1e4,
           mean = 10
  ),
  confparm('E_T',
           initial = 1,
           initial_condition_parameter = TRUE,
           prior = 'exponential',
           operator = 'realrw',
           estimate = TRUE,
           lower = 0,
           upper = 1e2,
           mean = 10
  ),
  confparm('T',
           initial = 1,
           initial_condition_parameter = TRUE,
           prior = 'exponential',
           operator = 'realrw',
           estimate = TRUE,
           lower = 0,
           upper = 1e2,
           mean = 10
  ),
  confparm('R_T',
           initial = 0,
           initial_condition_parameter = TRUE,
           prior = 'exponential',
           operator = 'realrw',
           estimate = TRUE,
           lower = 0,
           upper = 1e3,
           mean = 10
  ),
  confparm('S_H',
           initial = 250,
           initial_condition_parameter = TRUE,
           prior = 'exponential',
           operator = 'realrw',
           estimate = TRUE,
           lower = 0,
           upper = 1e3,
           mean = 10
  ),
  confparm('E_H',
           initial = 1,
           initial_condition_parameter = TRUE,
           prior = 'exponential',
           operator = 'realrw',
           estimate = TRUE,
           lower = 0,
           upper = 1e3,
           mean = 10
  ),
  confparm('H',
           initial = 1,
           initial_condition_parameter = TRUE,
           prior = 'exponential',
           operator = 'realrw',
           estimate = TRUE,
           lower = 0,
           upper = 1e3,
           mean = 10
  ),
  confparm('R_H',
           initial = 0,
           initial_condition_parameter = TRUE,
           prior = 'exponential',
           operator = 'realrw',
           estimate = TRUE,
           lower = 0,
           upper = 1e3,
           mean = 10
  )
)
```
```r
model <- config_phydyn(
  system.file('extdata', 'ba.2.86_algnWu-Hu-1.1_qc0_beast_template.xml', package = 'phydynbeast'),
  saveto = 'yor_output.xml',
  t0 = confparm('t0', estimate = FALSE, initial = 2021.9),
  equations = eqns,
  parameters = parms,
  coalescent_approximation = 'PL1',
  integrationSteps = 100,
  #method='adams-bashforth',
  minP = 0.001,
  penaltyAgtY = 0,
  useStateName = TRUE,
  traj_log_file = 'simodel0-traj.tsv',
  traj_log_frequency = 10000
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
