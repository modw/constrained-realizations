## Wiener filtering and constrained realizations with the Messenger Method
This code performs two tasks:
1. Given signal covariance (diagonal in ell space), noise covariance (diagonal in pixel space) and a mask, perform realizations of a fluctuation field. This field added to the Wiener filter of noisy data will give you constrained realizations of the data.
2. Given noisy data (d = s + n), mask and covariance matrices like the ones above, Wiener filter the data.

### Main Files:

* [`mm_constrained_realizations.py`](https://github.com/modw/constrained-realizations/blob/master/mm_constrained_realizations.py): Contains the `ConstrainedRealizations` class. The inputs are resolution, weights, frequency cut off and mask (optional). It has `set` methods for the signal and noise covariance matrices as well as the cooling schedule and random field seed. After the set-up, use `wiener_filter_data` to get Wiener Filter from map. The output is in harmonic space and should be transformed for visualization purposes. The `gen_constrained_realization` method generates one random fluctuation field that, added to the `wiener_filter`, will give you one random realization constrained by the data.
* [`support_classes.py`](https://github.com/modw/constrained-realizations/blob/master/support_classes.py): As the same suggests, it provides support classes to assist the `ConstrainedRealizations` class.


**Messenger Method**:

[1] Elsner, Franz, and Benjamin D. Wandelt. "Efficient Wiener filtering without preconditioning." Astronomy & Astrophysics 549 (2013): A111.

[2] Kodi Ramanah, Doogesh, Guilhem Lavaux, and Benjamin D. Wandelt. "Wiener filter reloaded: fast signal reconstruction without preconditioning." Monthly Notices of the Royal Astronomical Society 468.2 (2017): 1782-1793.
