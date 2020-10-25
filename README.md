# Megafauna Susceptibility

<strong>AUTHOR</strong>: Corey Bradshaw

<strong>CONTACT</strong>: corey.bradshaw@flinders.edu.au

<strong>URL</strong>: http://GlobalEcologyFlinders.com

<strong>INSTITUTION</strong>: Flinders University

<strong>RELEASE DATE</strong>: August 2020

R code accompanies article currently in review: 

BRADSHAW, CJA, CN JOHNSON, J LLEWELYN, V WEISBECKER, G STRONA, F SALTRÉ. In review. Relative demographic susceptibility does not explain the extinction chronology of Sahul’s megafauna. <i>eLife</i> (currently available as a pre-print on <i>bioRχiv</i>: DOI:<a href="http://doi.org/10.1101/2020.10.16.342303">10.1101/2020.10.16.342303</a>)

<strong>AIM</strong>: construct plausible stochastic demographic models for main Sahul megafauna to determine relative demographic susceptibility to environmental change & novel predation (human) sources

Builds models for following groups/genera:
- VOMBATIFORM HERBIVORES: <i>Diprotodon</i> (†), <i>Palorchestes</i> (†), <i>Zygomaturus</i> (†), <i>Phascolonus</i> (†), <i>Vombatus ursinus</i>
- MACROPODIFORM HERBIVORES: <i>Protemnodon</i> (†), <i>Osphranter rufus</i>, <i>Sthenurus</i> (†), <i>Simosthenurus</i> (†), <i>Procoptodon</i> (†), <i>Metasthenurus</i> (†), <i>Notamacropus</i>
- LARGE BIRDS: <i>Genyornis</i> (†), <i>Dromaius novaehollandiae</i>, <i>Alectura lathami</i>
- CARNIVORES: <i>Sarcophilus</i>, <i>Thylacinus</i> (†), <i>Thylacoleo</i> (†), <i>Dasyurus</i>
- MONOTREMES: <i>Megalibgwilia</i> (†), <i>Tachyglossus</i>

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
