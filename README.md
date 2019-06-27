# Predictive-Control
Some scripts connected with predictive control of two-dimensional neutralisation reactor. Scripts were written
for my master thesis:
## NONLINEAR MODEL PREDICTIVE CONTROL OF MULTIVARIABLE PROCESS USING A NEURAL WIENER STRUCTURE
Scripts covers:
* reactor modelling,
* creating data sets for training and validation purposes,
* linear and nonlinear models training
* predictive control of the model using MPC-NO, MPC-NPSL and GPC algorithms
* algorithms comparison

There is no master file to rule them all so scritps have to be used in the following way:
1. generate data sets (already generated - check "Dane" folder)
2. train models (already trained - check "Dane" folder)
3. run one of the "run*.m" for choosen algorithm

# Sources
* _P. Okoński (2018) Nonlinear model predictive control of multivariable process using a neural Wiener structure._

The thesis paper can be found in the Warsaw's University of Technology Base of Knowledge.

Description of the process can also be found in the following paper:
* _M. Ławryńczuk (2010) Suboptimal nonlinear predictive control based on multivariable neural Hammerstein models. M. Appl Intell 32:173-192_
