# Megafauna Susceptibility

<strong>AUTHOR</strong>: Corey Bradshaw

<strong>CONTACT</strong>: corey.bradshaw@flinders.edu.au

<strong>URL</strong>: http://GlobalEcologyFlinders.com

<strong>INSTITUTION</strong>: Flinders University

<strong>RELEASE DATE</strong>: August 2020

R code accompanies article currently in review: 

<a href="http://www.flinders.edu.au/people/corey.bradshaw">BRADSHAW, CJA</a>, <a href="https://www.utas.edu.au/profiles/staff/biological-sciences/chris-johnson">CN JOHNSON</a>, <a href="http://www.flinders.edu.au/people/john.llewelyn">J LLEWELYN</a>, <a href="https://researchnow.flinders.edu.au/en/persons/vera-weisbecker">V WEISBECKER</a>, <a href="https://researchportal.helsinki.fi/en/persons/giovanni-strona">G STRONA</a>, <a href="http://www.flinders.edu.au/people/frederik.saltre">F SALTRÉ</a>. In review. Relative demographic susceptibility does not explain the extinction chronology of Sahul’s megafauna. <i>eLife</i> (currently available as a pre-print on <i>bioRχiv</i>: DOI:<a href="http://doi.org/10.1101/2020.10.16.342303">10.1101/2020.10.16.342303</a>)

<strong>AIM</strong>: construct plausible stochastic demographic models for main Sahul megafauna to determine relative demographic susceptibility to environmental change & novel predation (human) sources

Builds models for following groups/genera:
- VOMBATIFORM HERBIVORES: <i>Diprotodon</i> (†), <i>Palorchestes</i> (†), <i>Zygomaturus</i> (†), <i>Phascolonus</i> (†), <i>Vombatus ursinus</i>
- MACROPODIFORM HERBIVORES: <i>Protemnodon</i> (†), <i>Osphranter rufus</i>, <i>Sthenurus</i> (†), <i>Simosthenurus</i> (†), <i>Procoptodon</i> (†), <i>Metasthenurus</i> (†), <i>Notamacropus</i>
- LARGE BIRDS: <i>Genyornis</i> (†), <i>Dromaius novaehollandiae</i>, <i>Alectura lathami</i>
- CARNIVORES: <i>Sarcophilus</i>, <i>Thylacinus</i> (†), <i>Thylacoleo</i> (†), <i>Dasyurus</i>
- MONOTREMES: <i>Megalibgwilia</i> (†), <i>Tachyglossus</i>

Repository includes the following files:

- '<a href="https://github.com/cjabradshaw/MegafaunaSusceptibility/blob/master/Sahul%20megafauna%20demographic%20susceptibility-base%20models.R">Sahul megafauna demographic susceptibility-base models.R</a>' — constructs bases models for all perturbation scenarios (must be run first)
- '<a href="https://github.com/cjabradshaw/MegafaunaSusceptibility/blob/master/matrixOperators.r">matrixOperators.R</a>' — functions to manipulate matrix models
- '<a href="https://github.com/cjabradshaw/MegafaunaSusceptibility/blob/master/megsuscept.SCENARIO2.juvsurv.R">megsuscept.SCENARIO2.juvsurv.R</a>' — runs Scenario 2 (reduction in juvenile survival)
- '<a href="https://github.com/cjabradshaw/MegafaunaSusceptibility/blob/master/megsuscept.SCENARIO3.fertred.R">megsuscept.SCENARIO3.fertred.R</a>' — runs Scenario 3 (reduction in fertility)
- '<a href="https://github.com/cjabradshaw/MegafaunaSusceptibility/blob/master/megsuscept.SCENARIO4.survred.R">megsuscept.SCENARIO4.survred.R</a>' — runs Scenario 4 (reduction in all-ages survival)
- '<a href="https://github.com/cjabradshaw/MegafaunaSusceptibility/blob/master/megsuscept.SCENARIO5.indrem.R">megsuscept.SCENARIO5.indrem.R</a>' — runs Scenario 5 (increasing individual offtake)
- '<a href="https://github.com/cjabradshaw/MegafaunaSusceptibility/blob/master/megsuscept.SCENARIO6.catinc.R">megsuscept.SCENARIO6.catinc.R</a>' — runs Scenario 6 (increasing frequency of catastrophic die-offs)
- '<a href="https://github.com/cjabradshaw/MegafaunaSusceptibility/blob/master/megsuscept.SCENARIO7.catMinc.R">megsuscept.SCENARIO7.catMinc.R</a>' — runs Scenario 7 (increasing magnitude of catastrophic die-offs)

** NOTE: For Scenario 7, Scenario 6 must be run first to create input .csv file 'catincQpr.csv' **

- '<a href="https://github.com/cjabradshaw/MegafaunaSusceptibility/blob/master/GRIWM.jk.sensitivity.R">GRIWM.jk.sensitivity.R</a>' — is the Gaussian-Resampled, Inverse-Weighted McInerny (GRIWM) Signor-Lipps algorithm, including a jack-knife estimator to test senstivity of dates to different assumptions. This R code is applied to the various chronologies for the following taxa (.csv files in the 'chronologies' sub-folder): <i>Diprotodon</i>, <i>Palorchestes</i>, <i>Zygomaturus</i>, <i>Phascolonus</i>, <i>Protemnodon</i>, <i>Sthenurus</i>, <i>Simosthenurus</i>, <i>Procoptodon</i>, <i>Metasthenurus</i>, <i>Genyornis</i>, <i>Thylacoleo</i>, <i>Megalibgwilia</i>
