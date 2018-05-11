# FBIM
Fluctuating Boundary Integral Method (FBIM)  by Yuanxun Bao, Manas Rachh, Eric E. Keaveny, Leslie Greengard and Aleksandar Donev

This repository contains demo codes for simulating the overdamped Brownian Dynamics of suspension of rigid particles using  boundary integral method. For details of the method, see our paper: [https://arxiv.org/abs/1709.01480](https://arxiv.org/abs/1709.01480)

## Table of Contents
libFBEM: contains routines needed for FBIM

libMatCode: contains routines for exporting matrices needed for FBIM

SingQuad: codes for exporting singular quadrature

SingleBodyTest: demo codes for BD simulation of a single body (disk, ellipse, starfish)

## Instructions
1. Open MATLAB and compile mex codes:
..* `cd libFBEM` and `mex -v fastgridding2d_mex.c`
..* `cd libMatCode` and `mex -v expint_eone.C -I/usr/local/include/ -lgsl`

2. 
