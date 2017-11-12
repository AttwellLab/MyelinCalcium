# MyelinCalcium
Code used for analysis of data in Krasnow, Ford, Valdivia, Wilson &amp; Attwell (2017) Nature Neuroscience

Run Code_transient_and_parameters.m (in MATLAB files folder) to plot calcium imaging data from GECIquant results for ROIs in sheaths, processes and somata.

Code requires the following functions:

Intersections
Version: 1.12, 27 January 2010
Author:  Douglas M. Schwarz

function [ycorr,yfit] = bf(y,varargin)

Functions by Anna Krasnow

read_txt.m

read_num.m

remove_spikes.m

Place xls files to be processed in the 'results files' folder.

Example file: GC1 Results, tab one with data from GeciQuant, tab two with false positives identified for removal. 

Use the same format for your files or update accordingly in the code.

Folder test_test contains results for when you run the code on the example file.

Adjust parameters in section '%% define variables'as required for your calcium traces.

Change date and folder names:

32/ date = '_test'; % change this date to match the analysed experiment

53/ mainFOLDER = strcat('Test', date);
