# Identity Exchange Profiler (IEP)

## How to install
```
pip install git+https://github.com/BartBruininks/IdentityExchangeProfiler
```

## How to use
Using segmentation data, like the leaflet identity per lipid over time (with or without noise) we can easily obtain the lipid flip-flop in an efficient and precise manner. We should specify which two labels we are interested in, in the case of the lipid flip-flop this would be the label for the outer and inner lealets. By default noise cancelling is used to plot the flip traces over time.

# Merging some compartments (optional)
```
selected_segmentation = np.copy(test_segmentation)

segment_A = 1 # inner leaflet
segment_B = 2 # outer leaflet
dt = 1000.0 # in ps, we save a frame every 1 ns

# Run analysis
analysis = iep.Profiler(
    segmentation=selected_segmentation, # selected_segmentation or test_segmentation
    segment_A=segment_A,
    segment_B=segment_B,
    dt = dt,
    verbose = True,
)
```
