% parameters seeting for calculating OBS tilt and compliance
% Modified by Zach Eilon 03/2016

clear all;
close all;

% path for matlab codes and functions
addpath ('function');

% location of the continures sacdata for event based
WORKINGdir = '/data/irma2/helenj/cascadia/EARTHQUAKES/SURFWAVES_DAYDATA/datacache_day/';

% EVTsacdata_input, SACPZ_input directory 
INSTRUMENTdir =  'NONE'; %'../INSTRUMENT/';
EVTsacdir     = '/data/irma2/helenj/cascadia/EARTHQUAKES/SURFWAVES_DATA/datacache/';


% channel naming
% chz='HZ';            %Channel name of Z component
% ch1='H1';    %Channel name of H1 component 
% ch2='H2';    %Channel name of H2 component
% chp='DH';	  %Channel name of DPG component

% day info to calculate tilt and compliance
nday = 2;     % how many days before event time.
T    = 6000;  % the legnth of each time window, in sec  
nwin = 20;    % numbers of time window in a day


fmin=0.005;fmax=0.02;
coh_min=0.8;
chorz_max=1;