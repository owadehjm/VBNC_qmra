# VBNC_qmra
This is a quantititative microbial risk assessment (QMRA) that accounts the formation of viable but nonculturalbe cells in Shiga toxin-producing Escherichia coli O157:H7 in lettuce during forward processing delays. The model conducts a farm-to-fork QMRA for STEC o157:H7 in lettuce for shredded and packaged lettuce (undergoes processing in  a plant) and field-bagged Romaine hearts at 10^5 Monte Carlo simulations.

bagged.code.R It is an R source code with the R function function called bagged_scenario that conducts a farm to form risk assessment for field-bagged Romaine hearts accounting for VBNC formation in cold storage during forward processing delays. Ann output of risk of illness per person per serving and annual risk of illness per person are generated as outputs. Sensitivity analysis is done using random forest function generating tornado plots.

computation.R. R script that runs the source codes for the field-bagged Romaine hearts and the shredded and packaged lettuce.

par_C.csv. A CSV file that provides the paired k1 and k2 parameters for the Juneja and Marks model for quantifying inactivation during chlorine washing process for the shredded and packaged lettuce. Find a separate analysis paper by Owade et al. 2024 detailing model performance in estimating STEC O157: H7 inactivation.

par_cool.csv. A CSV file that provides the paired k1 and k2 parameters for the Juneja and Marks model for quantifying inactivation at the cooling step during forward processing.

par_harvest.csv. A CSV file that provides the paired k1 and k2 parameters for the Juneja and Marks model for quantifying inactivation during harvesting of produce.

par_irrigation.csv. A CSV file that provides the paired k1 and k2 parameters for the Juneja and Marks model for quantifying inactivation during pre-harvsrvest holding time.

par_sto.csv. A CSV file that provides the paired k1 and k2 parameters for the Juneja and Marks model for quantifying inactivation during produce storage post processing and distribution.

processed.code.R. It is an R source code with the R function function called processed_scenario that conducts a farm to form risk assessment for shredded and packaged lettuce accounting for VBNC formation in cold storage during forward processing delays. An output of risk of illness per person per serving and annual risk of illness per person are generated as outputs. Sensitivity analysis is done using random forest function generating tornado plots.

Retail.temp.csv. A CSV file with empirical values with time-temperature profiles for post-processing and distribution of lettuce.
