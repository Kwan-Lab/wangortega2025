# Code For Analyzing Photostimulation Experiments

Frontal noradrenergic and cholinergic transients exhibit distinct spatiotemporal dynamics during competitive decision-making
Hongli Wang, Heather K. Ortega, Emma B. Kelly, Jonathan Indajang, Jiesi Feng, Yulong Li, Alex C. Kwan
doi: https://doi.org/10.1101/2024.01.23.576893

## System Requirements

- Tested on - Windows 10
-- 12th Gen Intel(R) Core(TM) i9-12900K, 3200 Mhz, 16 Core, 24 Logical Processor
-- 65 GB RAM
##
- MATLAB 2021b
-- Image processing toolbox (11.4)
-- Curve fitting toolbox (3.6)
-- Optimization toolbox (9.2)
-- Communications toolbox (7.6)
-- Statistics and machine learning toolbox (12.2)
##
- Import and save files from GitHub


This text you see here is *actually- written in Markdown! To get a feel
for Markdown's syntax, type some text into the left window and
watch the results in the right.

## Installation Guide

Dillinger uses a number of open source projects to work properly:

- Install MATLAB and required add-ons
##
- Data must be stored in the following structure
-- .....\level1\region\data
-- example '...\1Region\LM2\data'
-- 1Region must include folders for all regions to be analyzed
##


## Demo

To create Fig 6 F&G
- open 'master_MP_STIM_combine_animals.m'
- change root_path (line 11) to one level above folders containing powers to be analyzed
- run (will take ~20 minutes)

To create Fig 6 I&J 
- open 'master_MP_STIM_combine_animals_1reg.m'
- change root_path (line 11) to be one level above folders containing regions to be analyzed
- run (will take ~10 minutes)

To create Fig 7 H-I
- open 'master_preference.m'
- change root_path (line 7) to one level above powers to be analyzed on the simple choice task
- run (will take ~1 minute)

To create Fig 7 D-G
- load session of interest
- run 'plot_preference_single_session(trialData)'
- will take < 1 minute
