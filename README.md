# PRF model
This repository creates stimuli and BOLD time series to test PRF analytical methods.  It will serve as a place where we create these stimuli based on known PRF models and see how well we can recover the PRF parameters from the synthetic test stimuli.

We expect to produce increasing complex simulations of the bold time series, adding different PRF model shapes, noise properties, eye movement properties, and HRF models.  We will use these time series of stimuli and BOLD signals as inputs to test different pRF analytical software models.  We hope others can take advantage of the repository and the data - which we will share on a Flywheel or other open site - to test their implementations as well.  

## Stimulus

In which we define for people how we write the stimulus

We start with existing stimuli from experiments we have run.

Then, we write code to create new stimuli, saving (a) the aperture with contrast stimuli carriers and (b) just the aperture (binary version with 1 for contrast and 0 elsewhere). This will let us explore the impact of bars of different widths, on/off-times, and other stimulus selection parameters. 

We also expect to implement rings and wedges, and we expect the KK style stimuli.

We will probably create these stimuli using some combination of vistadisp and the KK methods, but hopefully we can do most of it in vistadisp and stay with the vistasoft programmatic effort.  KK has many special and redundant routines.

## BOLD Time series

To write out the BOLD time series properties, we need to make some choices.

We need to choose the HRF model.  We need to choose the TR so the impact of stimulus and MR acquisition timing makes sense.

We also need to choose a noise model.  We might use Gaussian, or we might use Poisson.

We might incorporate the possibility that the pRF is jittered slightly over time, accounting for eye movements.

## PRF model

We start with estimates of the linear, circular PRF model.  

Then we will move to the CSS circular model (adds the exponent).  And then to the CSS with a potentially elliptical shape model.  

Finally, we expect to go all the way to SOC models that account for the carrier (contrast), not just the aperture.

We will evaluate the models by creating repeats of the same model with different noise samples.  We fit  one set of data with two different models, and we compare on a new batch of synthetic data with different noise samples. This is the usual cross-validation approach.

# Data sharing

We aren't yet sure where we will post these data and the algorithms.

Probably in Flywheel as data and Gears.

We may also put them on the Stanford Digital repository and Dockerhub.


