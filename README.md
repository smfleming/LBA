# LBA
Toolbox for fitting the Linear Ballistic Accumulator (LBA) model of decision-making (Brown & Heathcote, 2008, Cogn Psychology) in Matlab.
The code is adapated from Scott Brown's R code for fitting the LBA available here: http://www.newcl.org/brown

The code requires the Optimization Toolbox as maximum likelihood estimation relies on fmincon.

To get started, see LBA_example and LBA_wrapper.
Type "help LBA_mle" to see the full range of model fitting options. Through the "model" input the user can compare different models, e.g. by allowing only the drift rate or only the bound to vary across conditions. LBA_wrapper has an example of how to do this.

Please get in touch with any bugs, comments or questions at sf102@nyu.edu

