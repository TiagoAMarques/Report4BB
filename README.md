# Report4BB

This git repository has evolved over time and it now contains a set of analysis and explorations revolving arount the topic of how to estimate the precision on a prediction for a Gamma GLMM for unobserved levels of random effcts(s).

This endeavour originated a question on Stack Exchange, hosted here:

https://stats.stackexchange.com/questions/616697/how-to-estimate-precision-on-a-prediction-for-a-glmm-for-an-unobserved-level-of/

Said question let to some interaction with Ben Bolker (BB), who replied to it, and this repository was created first to share things with him.

After that, trying to implement the suggestions from BB, led me to a quite exaustive exploration of the differences between implementing the above model in lme4::glmer versus in a Bayesian context, in this case implemented in Nimble. 

This led to another question in the nimble users google group:

https://groups.google.com/g/nimble-users/c/7UHlKdCC8B4

The result of that is inside folder LookingAtNimble, and I try to reconcile all I found along the way in file "SpermWhaleCueRatesSuppPub1.Rmd".

Aknowledgements:

Ben Bolker, Ben Augustine, Perry de Valpine, Luca, Wei Zhang and Len Thomas provided feedback which might be reflected in this material. All pending errors are solely mine. 
