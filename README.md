# Megafauna Susceptibility

AUTHOR: Corey Bradshaw
CONTACT: corey.bradshaw@flinders.edu.au
URL: http://GlobalEcologyFlinders.com
INSTITUTION: Flinders University
RELEASE DATE: August 2020

R code accompanies article currently in review: 
BRADSHAW, CJA, CN JOHNSON, J LLEWELYN, V WEISBECKER, G STRONA, F SALTRÉ. In review. Relative demographic susceptibility does not explain the extinction chronology of Sahul’s megafauna.

AIM: construct plausible stochastic demographic models for main Sahul megafauna to determine relative demographic susceptibility to environmental change & novel predation (human) sources

Builds models for following groups/genera:
- VOMBATIFORM HERBIVORES: ✓Diprotodon (†), ✓Palorchestes (†), ✓Zygomaturus (†), ✓Phascolonus (†), ✓Vombatus ursinus
- MACROPODIFORM HERBIVORES: ✓Protemnodon (†), ✓Osphranter rufus, ✓Sthenurus (†), ✓Simosthenurus (†), ✓Procoptodon (†), ✓Metasthenurus (†), ✓Notamacropus
- LARGE BIRDS: ✓Genyornis (†), ✓Dromaius novaehollandiae, ✓Alectura lathami
- CARNIVORES: ✓Sarcophilus, ✓Thylacinus (†), ✓Thylacoleo (†), ✓Dasyurus
- MONOTREMES: ✓Megalibgwilia (†), ✓Tachyglossus

Repository includes the following files:

- 'Sahul megafauna demographic susceptibility-base models.R' — constructs bases models for all perturbation scenarios (must be run first)
- 'matrixOperators.R' — functions to manipulate matrix models
- 'megsuscept.SCENARIO2.juvsurv.R' — runs Scenario 2 (reduction in juvenile survival)
- 'megsuscept.SCENARIO3.fertred.R' - runs Scenario 3 (reduction in fertility)
- 'megsuscept.SCENARIO4.survred.R' - runs Scenario 4 (reduction in all-ages survival)
- 'megsuscept.SCENARIO5.indrem.R' - runs Scenario 5 (increasing individual offtake)
- 'megsuscept.SCENARIO6.catinc.R' - runs Scenario 6 (increasing frequency of catastrophic die-offs)
- 'megsuscept.SCENARIO7.catMinc.R' - runs Scenario 7 (increasing magnitude of catastrophic die-offs)

** NOTE: For Scenario 7, Scenario 6 must be run first to create input .csv file 'catincQpr.csv' **
