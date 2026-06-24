# Catenary Curve: Variational Derivation and Newton Solver

This repository contains a small educational Python project devoted to the calculation of a catenary curve suspended between two fixed points. The problem is treated as a classical example of the variational method: the equilibrium shape of the curve is obtained by minimizing its gravitational potential energy under the constraint of fixed length.

Although the problem is mathematically simple, it is useful from both pedagogical and practical points of view. It demonstrates how a physical equilibrium problem can be formulated as a constrained variational problem, reduced to a nonlinear equation, and then solved numerically.

Educational purpose

The project illustrates several important ideas:

formulation of the catenary problem as minimization of a functional;
use of a Lagrange multiplier to account for the fixed-length constraint;
derivation of the catenary equation from the Euler equation;
reduction of the boundary-value problem to a single nonlinear equation;
numerical solution of this equation by Newton’s method;
calculation of the tension force and tensile stress along the curve;
visualization of the catenary shape and stress distribution using Python.
Physical problem

A flexible curve with uniform linear mass density is suspended by its fixed endpoints in a gravitational field. The total length of the curve is prescribed and must be larger than the straight-line distance between the endpoints.

The goal is to determine:

the equilibrium shape of the curve;
the horizontal and vertical components of the tension force;
the total tangential tension;
the tensile stress in the material.

The endpoints may be located at different heights.

Numerical method

After introducing suitable dimensionless variables, the problem is reduced to a nonlinear equation for a single unknown parameter. This equation is solved using Newton’s method. The algorithm also takes into account the removable singularity that appears in the limiting case of small sag.

Python implementation

The Python script calculates:

parameters of the catenary curve;
minimum admissible length;
curve shape;
tension force along the curve;
tensile stress along the curve;
minimum and maximum values of tension and stress.

The script also produces plots of:

the calculated catenary curve with fixed endpoints;
the tension and tensile stress along the curve.
Practical relevance

The maximum tensile stress is an important engineering parameter when selecting the cross-section and length of a suspended wire, cable, or conductor. In practical applications, it should be compared not directly with the yield strength of the material, but with the allowable tensile stress defined using a safety factor.

Requirements

The script requires:

numpy
matplotlib
Running the script

Run the Python file:

python catenary_curve.py

The program asks for the total curve length and then calculates the catenary parameters, tension, stress, and plots the results.

Repository status

This is an educational project. It is intended for learning, demonstration, and simple engineering estimates rather than for certified structural or mechanical design.
