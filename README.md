# DistributedFastTopKSelection
Matlab Code for "Fast networked data selection via distributed smoothed quantile estimation"

## Distributed smoothed quantile estimation

#### Algorithm comparison for different cases
- Main.m: normal case, single data point within each node
- Main_LargeScale.m: large scale case, single data point within each node
- Main_MultNum.m: multiple data points within each node

#### Choice of smoothing techniques
- smooth='Nesterov';
- smooth='Convolution';

#### Choice of loss
- loss='l2';
- loss='l1';
- loss='inf';

#### ALgorithms
- DistributedQuantileEstimation_EXTRA.m: Distributed top-k selection via EXTRA for single data point within each node
- DistributedQuantileEstimation_EXTRA_MultNum.m: Distributed top-k selection via EXTRA for multiple data points within each node
- DistributedQuantileEstimation_SGD.m: Distributed top-k selection via DGD for single data point within each node
- DistributedQuantileEstimation_SGD_MultNum.m: Distributed top-k selection via DGD for multiple data points within each node
