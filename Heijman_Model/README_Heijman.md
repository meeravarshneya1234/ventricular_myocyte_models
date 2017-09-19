The original code for this model comes from the Rudy Lab website: http://rudylab.wustl.edu/research/cell/code/AllCodes.html.

Additions:
1) The runHeijman.m is a script file to make it easier to make revisions to settings and run the model. 
2) The determineAPD function within the mainHRdBA.m function was altered to replicate the find_APD.m function included in the main folder. 
It calculates APD as the as difference between time of stimulus and time when voltage returns to -75 mV. 
