clear all
close all
clc

format long
TDPR_data = load('TDPR.mat');
TWS_data = load('TWS.mat');

x = TDPR_data.x;

TDPR_P = TDPR_data.P;
TWS_P = TDPR_data.P;

(TDPR_P - TWS_P)'
