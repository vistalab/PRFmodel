# PRF model
This repository creates stimuli and BOLD time series to test PRF analytical methods.  It will serve as a place where we create these stimuli based on known PRF models and see how well we can recover the PRF parameters from the synthetic test stimuli.

We expect to produce increasing complex simulations of the bold time series, adding different PRF model shapes, noise properties, eye movement properties, and HRF models.  We will use these time series of stimuli and BOLD signals as inputs to test different pRF analytical software models.  We hope others can take advantage of the repository and the data - which we will share on a Flywheel or other open site - to test their implementations as well.  

## Stimulus

In which we define for people how we write the stimulus

We start with some existing stimuli.

Then, we create new stimuli (aperture with contrast stimuli carriers) and just the aperture (binary version with 1 for contrast and 0 elsewhere). This will let us have bars of different widths, off-times, and other parameters. We also expect rings and wedges, and we expect the KK style of put the image here and there.

We will probably try to do it with vistadisp in the long run.

## BOLD Time series

In which we define for people how we write out the BOLD time series properties.

We need to choose the HRF model as part of this.  Also the TR so the timing makes sense.

We also need to choose a noise model.  We might use Gaussian, or we might use Poisson.

## PRF model

We start with estimates of the either linear, circular PRF model.  We will quickly move to the CSS circular model (adds the exponent).  And then to the CSS with a potentially elliptical shape model.

We will be able to evaluate the models by creating repeats of the same model with different noise samples.  So we fit on one set of data with two different models, and then we compare on a new batch of data with different noise samples which model fits better.  The usual cross-validation approach.

# Data sharing

We aren't yet sure where we will post these data.  Probably in Flywheel, but we may also put them on the Stanford Digital repository or somewhere else.
