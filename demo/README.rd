what this module can do is ....


- set up the initial condition for a sequence of histones.(init_genome)
  - histones can be characterized with the following points.
    - histone status. M, U or A

    - histone kpluc value, kminus value
      - what is kplus?
      	- kplus is the prob. of getting an unmethylated histone methylated when the neigbors of the histone is methylated.
      - what is kminus?
      	- kminus is the prob. of turns a methylated histone unmethylated OR turns acetylated histone turns unmethylated.
    - histone knuc value, kace value
      - what is knuc?
      	- knuc is the prob. of getting E0 position turned methylated when transcription does not happen.
      - what is kace?
        - kace is the prob. of getting an unmethylated histone to turn acetylated.


- set up the environment that you want to calculate/visualize the sequence of histones.(track_epigenetic_process)
  - time for tracking (in hour scale)
  - window is the locus that we are interested in.

