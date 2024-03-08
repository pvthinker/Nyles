# Nyles

Nyles -- New Python LES (Large Eddy Simulation) -- is a code to simulate
incompressible flows.  It has several original features, especially its
way of handling the SGS (subgrid-scale) physics.  This is done by
upwinding the vorticity in the vortex force term, which shows up in the
covariant form of the momentum equation (Navier--Stokes equation).

The project started in September 2019 as a GFD project with two
graduate students of Marine Physics. We have now reached a point where
the numerical choices reveal their promises.

In 2023 WENO reconstructions was added, a new multigrid solver was added, and an attempt to implement the diagnostic of implicit dissipation.

To run it:

  - cd Nyles
  - make [for once, this compiles the Fortran into *.so libraries]
  - source activate.sh [every time you open a new terminal]
  - cd experiments/lock
  - python3 lockexchange.py

![ScreenShot](/screenshots/nyles_0.png)
