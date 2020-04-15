# Overview
Thus far, most studies for metabolic flux estimation have been done under steady state conditions. Few studies have been done under in vivo-like non-steady state conditions because of technical difficulties in both experiments and mathematical modeling. Furthermore, identification of regulation of glucose metabolism by phosphorylation events and allosteric effects is presently inferred from kinetics of isolated enzymes; consequently, which phosphorylation events and allosteric effectors regulate glucose flux in cells remains unknown. We solved these issues by developing novel mathematical approaches using data from phosphoproteomic and 13C-glucose-labelled metabolomic experiments. We applied the developed approach to insulin-stimulated adipocytes (Ohno et al., submitted)

Here we provide the MATLAB source code for 13C-metabolic flux analysis under non-steady-state conditions and kinetic modelling to identify key regulatory mechanisms for a reaction and relative contributions of these molecules to flux.

# Requirements
The developmental version of package has been tested on the following systems:

Windows 10

Excel 2016

MATLAB R2019a with Optimization Toolbox

# Set up
Download the package code to your local path (< 1 min)

# Contents
## MFA_NSS
This is the package code for metabolic flux analysis under non-steady-state.
Please run demo_MFA_NSS1.m and demo_MFA_NSS2.m in the package. Metabolic flux analysis for simple metabolic network will be performed. Your can also change population sizes and maximum nuber of function evaluations for metaheuristic optimization so that you can obtain better fitting results.

## KineticModeling
This is the package code for kinetic modeling of flux with regulatory molecules.
Please run demo_KM.m in the package. Kinetic modeling for phosphofructokinase 1 will be performed. Note that only a part of phosphorylation and allosteric effectors are considered as candidate of regulatory molecule for the reaction.

# Contact
Satoshi Ohno: sohno@bs.s.u-tokyo.ac.jp

Shinya Kuroda: skuroda@bs.s.u-tokyo.ac.jp

# Reference
Satoshi Ohno, Lake-Ee Quek, James R. Krycer, Katsuyuki Yugi, Akiyoshi Hirayama, Satsuki Ikeda, Futaba Shoji, Kumi Suzuki, Tomoyoshi Soga, David E. James, and Shinya Kuroda. Kinetic trans-omic analysis reveals key regulatory mechanisms for insulin-regulated glucose metabolism in adipocytes. Submitted.

