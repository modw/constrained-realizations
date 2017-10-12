## Wiener filtering and constrained realizations with the Messenger Method
This code performs two tasks:
1. Given signal covariance (diagonal in ell space), noise covariance (diagonal in pixel space) and a mask, perform realizations of a fluctuation field. This field added to the Wiener filter of noisy data will give you constrained realizations of the data.
2. Given noisy data (d = s + n), mask and covariance matrices like the ones above, Wiener filter the data.


**Messenger Method**:

[1] Elsner, Franz, and Benjamin D. Wandelt. "Efficient Wiener filtering without preconditioning." Astronomy & Astrophysics 549 (2013): A111.

[2] Kodi Ramanah, Doogesh, Guilhem Lavaux, and Benjamin D. Wandelt. "Wiener filter reloaded: fast signal reconstruction without preconditioning." Monthly Notices of the Royal Astronomical Society 468.2 (2017): 1782-1793.
