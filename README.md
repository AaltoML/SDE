# Applied Stochastic Differential Equations

[Simo Särkkä](https://users.aalto.fi/~ssarkka/) · [Arno Solin](http://arno.solin.fi)

Example codes for the book:

* Simo Särkkä and Arno Solin (2019). **Applied Stochastic Differential Equations**. Cambridge University Press. Cambridge, UK.

The book can be ordered through [Cambridge University Press](https://www.cambridge.org/fi/academic/subjects/statistics-probability/applied-probability-and-stochastic-networks/applied-stochastic-differential-equations?format=PB) or, e.g., from  [Amazon](https://www.amazon.co.uk/Stochastic-Differential-Equations-Mathematical-Statistics/dp/1316649466).

*With permission from the publisher, we are providing a PDF version of the book [here](https://users.aalto.fi/~asolin/sde-book/sde-book.pdf). This PDF version is made available for personal use. The copyright in all material rests with the authors (Simo Särkkä and Arno Solin). Commercial reproduction is prohibited, except as authorised by the author and publisher.*

## Summary

The book *Applied Stochastic Differential Equations* gives a gentle introduction to stochastic differential equations (SDEs). The low learning curve only assumes prior knowledge of ordinary differential equations and basic concepts of statistic, together with understanding of linear algebra, vector calculus, and Bayesian inference. The book is mainly intended for advanced undergraduate and graduate students in applied mathematics, signal processing, control engineering, statistics, and computer science.

The worked examples and numerical simulation studies in each chapter illustrate how the theory works in practice and can be implemented for solving the problems. To promote hands-on work with the methods, we provide the MATLAB􏰀 source code for reproducing the example results in the book. The code examples have been grouped by chapter, and some pointers to example and figure numbers in the book are given below.

## Codes for the examples

All codes for the examples (excluding pen and paper examples which do not have any code associated to them) and figures (those requiring numerical simulation) in the book are provided below. Some entries cover multiple examples/figures in the same chapter, in the case of which the code file naming follows the first example.

### Chapter 2: Numerical solution of ODEs

This experiment replicates the results in Example 2.9  (Numerical solution of ODEs) in the book (Fig. 2.1).

```ch02_ex09_numerical_solution_of_odes.m```

### Chapter 3: Two views of Brownian motion 
  
This experiment replicates the results in Figure 3.2 in the book.

```ch03_fig02_two_views_of_brownian_motion.m```
  
### Chapter 3: Time-varying oscillator 

This experiment replicates the results in Example 3.7  (Heart and breathing tracking in the brain) in the book (Fig. 3.8).

```ch03_ex07_time_varying_oscillator.m```

### Chapter 3: Stochastic spring model
  
This experiment replicates the results in Example 3.10 (Stochastic spring model) in the book (Figs. 3.9 and 3.11).

```ch03_ex10_stochastic_spring_model.m```
  
### Chapter 3: White noise

This experiment replicates the results in Figure 3.10 in the book.

```ch03_fig10_white_noise.m``` 

### Chapter 4: Brownian motion

This experiment replicates the results in Figure 4.1 in the book.

```ch04_fig01_brownian_motion.m```
  
### Chapter 4: Solution of the Ornstein–Uhlenbeck process

This experiment replicates the results in Example 4.5 (Solution of the Ornstein–Uhlenbeck process) in the book (Fig. 4.2).

```ch04_ex05_ou_process.m```
  
### Chapter 7: Karhunen–Loeve series of Brownian motion

This experiment replicates the results in Example 7.3 (Karhunen–Loeve series of Brownian motion) in the book (Fig. 7.1).

```ch07_ex03_karhunen_loeve_series.m```

### Chapter 7: Doob's h-transform  
  
This experiment replicates the results in Example 7.12 (Conditioned Ornstein–Uhlenbeck process) in the book (Fig. 7.2).

```ch07_ex12_doobs_h_transform.m```    

### Chapter 7: Feynman-Kac formula

This experiment replicates the results in Example 7.17 (Solution of an elliptic PDE using SDE simulation) in the book (Fig. 7.3).

```ch07_ex17_feynman_kac.m```
  
### Chapter 8: Weak Gaussian vs. weak three-point approximations  

This experiment replicates the results in Example 8.6 (Simulating from a trigonometric nonlinear SDE) in the book (Fig. 8.1).

```ch08_ex06_weak_itotaylor.m```
  
### Chapter 8: Comparison of Runge–Kutta schemes

This experiment replicates the results in Example 8.11 (Comparison of ODE solvers) in the book (Fig. 8.2).

```ch08_ex11_ode_solvers.m``` 

### Chapter 8: Duffing van der Pol model

This experiment replicates the results in Example 8.15 (Duffing van der Pol oscillator) in the book (Figs. 8.3–8.5).

```ch08_ex15_duffing_van_der_pol.m```
  
### Chapter 8: Leapfrog Verlet

This experiment replicates the results in Example 8.21 (Leapfrog solution to the spring model) in the book (Fig. 8.6).

```ch08_ex21_leapfrog_verlet.m```    
  
### Chapter 8: Exact simulation of sine diffusion

This experiment replicates the results in Example 8.24 (Exact simulation of sine diffusion) in the book (Fig. 8.7).

```ch08_ex24_exact_simulation.m```    
  
### Chapter 9: Linearizations and approximations for the Beneš model

This experiment replicates the results in Example 9.6 (Linearization and Gauss–Hermite approximations, shown in Fig. 9.1) and Example 9.14 (Local linearization vs. Gaussian approximations, shown in Fig. 9.3).

```ch09_ex06_linearizations_for_benes.m```
  
### Chapter 9: Gaussian approximation

This experiment replicates the results in Example 9.7 (Gaussian approximation of a nonlinear trigonometric SDE) in the book (Fig. 9.2).

```ch09_ex07_gaussian_approximation.m```
  
### Chapter 9: Hermite expansion of Beneš SDE 

This experiment replicates the results in Example 9.18 (Hermite expansion of Beneš SDE) in the book (Fig. 9.4).

```ch09_ex18_hermite_expansion.m```
  
### Chapter 9: Discretized FPK for the Beneš SDE  

This experiment replicates the results in Example 9.20 (Discretized FPK for the Beneš SDE) in the book (Fig. 9.5).

```ch09_ex20_discretized_fpk_for_benes.m```
  
### Chapter 9: Pathwise series expansion of the Beneš SDE

This experiment replicates the results in Example 9.23 (Pathwise series expansion of Beneš SDE) in the book (Fig. 9.6).

```ch09_ex23_pathwise_series_expansion.m```

### Chapter 10: Beneš-Daum and EKF/ERTS examples

This experiment replicates the results in Example 10.17 (Beneš–Daum filter, Fig. 10.2), Example 10.26 (Continuous-discrete EKF solution to the Beneš–Daum problem, Fig. 10.4), and  Example 10.38 (Extended RTS solution to Beneš and Beneš–Daum filtering problems, Fig. 10.6) in the book.

```ch10_ex17_benes_daum_ekf_erts.m```  
  
### Chapter 10: Ornstein-Uhlenbeck filtering and smoothing

This experiment replicates the results in Example 10.19 (Kalman filter for the Ornstein–Uhlenbeck model) and Example 10.21 (Continuous-discrete Kalman filter for the Ornstein–Uhlenbeck model), the results of which are shared in Figure 10.3. Additionally, the experiment also covers Example 10.29 (RTS smoother for the Ornstein–Uhlenbeck model) and Example 10.33 (Continuous-discrete/continuous RTS smoother for the Ornstein–Uhlenbeck model), the results of which are in Figure 10.5.

```ch10_ex19_ou_filtering_smoothing.m```  
  
### Chapter 11: Parameter estimation in Ornstein-Uhlenbeck

This experiment replicates the results in Example 11.5 (Exact parameter estimation in the Ornstein–Uhlenbeck model, Fig. 11.1) and Example 11.9 (Approximate parameter estimation in the Ornstein–Uhlenbeck model, Fig. 11.2) in the book.

```ch11_ex05_ou_parameter_estimation.m```
  
### Chapter 12: Batch and sequential solution to GP regression
  
This experiment replicates the results in Example 12.6 (Batch GP regression, Fig. 12.1) and Example 12.11 (Sequential solution to GP regression, Fig. 12.2) in the book.

```ch12_ex06_gp_regression.m```  
  
### Chapter 12: GP approximation of the drift function (double-well)

This experiment replicates the results in Example 12.12 (GP approximation of the double-well model) in the book (Fig. 12.3).

```ch12_ex12_gp_approximation_of_drift.m```

## How to run?

These codes have been tested under *Mathworks MATLAB R2018b* and *GNU Octave 4.4*. The code is pseudo-code like and tries to follow the presentation in the book. Thus they should be applicable for porting to other languages as well (consider this an exercise).

## License

Copyright Simo Särkkä and Arno Solin.

This software is provided under the MIT License. See the accompanying LICENSE file for details.
