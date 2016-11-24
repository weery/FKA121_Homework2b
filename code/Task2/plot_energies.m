clear all, close all, clc

energies=importdata('energy.dat');

x=(1:length(energies))';
mean = cumsum(energies)./x;

plot(mean)