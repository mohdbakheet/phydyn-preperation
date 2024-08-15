# Using phydynBEAST for Multi-Host SIR Models

This README provides instructions on how to use the `phydynBEAST` package developed by Erik Volz to generate BEAST2 XML files for a Susceptible-Infected-Recovered (SIR) compartmental model across multiple host species.

## Installation

1. **Install R:** Make sure you have R installed on your computer (download from https://cran.r-project.org).

2. **Install phydynBEAST:**
   ```R
   devtools::install_github('emvolz/phydynbeast')
   ```

## Setting Up Your Model

1. **Creating xml file
   The first step is by generating an xml templete using a Yule tree prior and then removing xml elements that specify the Yule prior and operator.
3. **Load the package:**
   ```R
   library(phydynBEAST)
   ```

4. **Define the multi-host SIR model:**
   ```R
  # Define the states for the SIR model
   eqns <- list( 
  confeqn('N_C = S_C + E_C + C + R_C', type='definition'),
  confeqn('N_T = S_T + E_T + T + R_T', type='definition'),
  confeqn('N_H = S_H + E_H + H + R_H', type='definition'),
 
  
  confeqn('beta_CC * C / N_C', origin='C', type='birth', destination='E_C'),
  confeqn('beta_TC * T / N_T', origin='T', type='birth', destination='E_T'),
  confeqn('beta_HC * H / N_H', origin='H', type='birth', destination='E_H'),
  
  
  confeqn('beta_CT * C / N_C', origin='C', type='birth', destination='E_db'),
  confeqn('beta_TT * T / N_T', origin='T', type='birth', destination='E_db'),
  confeqn('beta_HT * H / N_H', origin='H', type='birth', destination='E_db'),
 
  confeqn('beta_CH * C / N_C', origin='C', type='birth', destination='E_wm'),
  confeqn('beta_TH * T / N_T', origin='T', type='birth', destination='E_wm'),
  confeqn('beta_HH * H / N_H', origin='H', type='birth', destination='E_wm'),
  
  confeqn('(gamma_C + mu_H) * C',  origin='C', type='death'),
  confeqn('(gamma_T + mu_H) * T',  origin='T', type='death'),
  confeqn('(gamma_H + mu_H) * H',  origin='H', type='death'),
  
  confeqn('eta_C * E_T', origin='E_C', type='migration', destination='C'),
  confeqn('eta_T * E_T', origin='E_T', type='migration', destination='T'),
  confeqn('eta_H * E_H', origin='E_H',  type='migration', destination='H'),
  
  confeqn('gamma_C * C', origin='R_C', type='nondeme'),
  confeqn('gamma_T * T', origin='R_T', type='nondeme'),
  confeqn('gamma_H * H', origin='R_H',  type='nondeme' ),

  confeqn('-(beta_CC * C / N_C + beta_TC * T / N_T + beta_HC * H / N_H) * S_C', origin='S_C', type='nondeme'),
  confeqn('-(beta_CT * C / N_C + beta_TT * T / N_T + beta_HT * H / N_H) * S_T', origin='S_db', type='nondeme'),
  confeqn('-(beta_CH * C / N_C + beta_TH * T / N_T + beta_HH * H / N_H) * S_H', origin='S_h', type='nondeme')
  
)
   params <- c("beta_wb", "gamma_wb", "beta_db", "gamma_db")

   # Define prameters
  

   # Create the model
   ```

5. **Set initial conditions and other parameters:**
   ```R
   # Define initial conditions
  
   ```

## Generating the BEAST2 XML File

1. **Generate XML:**
   ```R
   ```

2. **Save the XML to a file:**
   ```R
   ```

## Running the Analysis

After generating the XML file, you can run the analysis using BEAST2 software. Make sure you have BEAST2 installed (download from http://beast2.org).

1. **Run BEAST2:**
   - Open BEAST2 and select the `SIR_MultiHost.xml` file.
   - Configure your run settings.
   - Run the analysis.

## Conclusion

This setup provides a basic framework for modeling multi-host SIR dynamics using the `phydynBEAST` package. Adjust the model and parameters according to your specific research needs.
