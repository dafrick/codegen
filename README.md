# Automatic code-generation tool for hybrid MPC
a low-complexity, iterative method for hybrid model predictive control 
## Introduction
In [\[arXiv:1609.02819\]](http://arxiv.org/abs/1609.02819) we present a novel iterative method for finding loval minima of hybrid MPC roblems 
## Installation instructions
Make sure you have installed the [multi-parametric toolbox](http://control.ee.ethz.ch/~mpt/3/Main/Installation) (MPT3)
## How to use
### Example
A simple example (example.m) illustrates the use of the tool. It can be run by typing ```example()``` while in the directory of the tool.
### Generate solver
```pcg.generateSolver()``` is the function that generates a parametric solver for a given problem. Just type ```help pcg.generateSolver``` in MATLAB to find out more.
