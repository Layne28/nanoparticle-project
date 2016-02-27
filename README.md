The goal of this project is to simulate the chemical etching of a gold nanoparticle. The project was inspired by work done in Paul Alivisatos' group at Berkeley.
My approach, at least initially, has been to use a Grand Canonical Metropolis Monte Carlo (GCMMC) algorithm to simulate the etching process. Specifically, atom
insertion and deletion moves are generated and then either accepted or rejected according to the values of the specified chemical potential and bond energy.
