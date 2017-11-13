# MyelinCalcium
Code used for analysis of data in Krasnow, Ford, Valdivia, Wilson &amp; Attwell (2017) Nature Neuroscience

Run Code_transient_and_parameters.m (in MATLAB files folder) to plot calcium imaging data from GECIquant results for ROIs in sheaths, processes and somata.

This code requires the following functions:

(1) Intersections (MatLab, included in this repository)
Version: 1.12, 27 January 2010
Author:  Douglas M. Schwarz

(2) function [ycorr,yfit] = bf(y,varargin)  (also in this repository)

MatLab Functions by Anna Krasnow  (also in this repository)

read_txt.m

read_num.m

remove_spikes.m

What you should do:

Place your data (xls) files to be processed in the 'results files' folder.

An example file is provided: GC1 Results, tab one with data from GeciQuant, tab two with false positives identified for removal. 

Use the same format for your files or update accordingly in the code.

The folder test_test (with the format FolderName_Date) contains results for when you run the code on the example file.

Adjust the parameters in the section '%% define variables'as required for your calcium(t) traces.

Change the date and folder names:

32/ date = '_test'; % change this date to match the analysed experiment

53/ mainFOLDER = strcat('Test', date); % change this folder to a name for the analysed experiment
