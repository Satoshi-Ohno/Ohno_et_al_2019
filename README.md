# Overview
Thus far, most studies for metabolic flux estimation have been done under steady state conditions. Few studies have been done under in vivo-like non-steady state conditions because of technical difficulties in both experiments and mathematical modeling. Furthermore, identification of regulation of glucose metabolism by phosphorylation events and allosteric effects is presently inferred from kinetics of isolated enzymes; consequently, which phosphorylation events and allosteric effectors regulate glucose flux in cells remains unknown. We solved these issues by developing novel mathematical approaches using data from phosphoproteomic and 13C-glucose-labelled metabolomic experiments. We applied the developed approach to insulin-stimulated adipocytes (Ohno et al., 2019)

Here we provide the MATLAB source code for metabolic flux analysis under non-steady-state and identification of flux regulators which contribute to the change in flux through a particular metabolic reaction. Such regulators were denoted “effective flux regulators.”

# Requirements
The developmental version of package has been tested on the following systems:
Windows 10
MATLAB R2017a with Optimization Toolbox

# Set up
Download the package code to your local path (< 1 min)

# Contents
MetabolicFluxAnalysis_NonSteadyState
Package code for metabolic flux analysis under non-steady-state.
Please run demo_MFA_NSS.m in the package. Metabolic flux analysis for simple metabolic network will be performed.

IdentificaitonIOfEffectiveFluxRegulators
Package code for Identification of effective flux regulators.
Please run demo_IEFR.m in the package. Identification of effective flux regulators for glycogen synthase will be performed. Note that only a part of phosphorylation and allosteric effectors are considered as candidate of effective flux regulators.

# Contact
Shinya Kuroda: skuroda@bs.s.u-tokyo.ac.jp

# Reference
Satoshi Ohno, Lake-Ee Quek, James R. Krycer, Katsuyuki Yugi, Akiyoshi Hirayama, Satsuki Ikeda, Futaba Shoji, Kumi Suzuki, Tomoyoshi Soga, David E. James, and Shinya Kuroda. Key regulatory mechanisms for insulin-induced changes in glucose metabolism in adipocytes. Submitted.

