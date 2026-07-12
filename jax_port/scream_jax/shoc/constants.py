"""SHOC-specific constants.

Transcribed verbatim from
components/eamxx/src/physics/shoc/shoc_constants.hpp
(struct scream::shoc::Constants<Scalar>).
"""

mintke = 0.0004          # Minimum TKE [m2/s2]
maxtke = 50.0            # Maximum TKE [m2/s2]
minlen = 20.0            # Lower limit for mixing length [m]
maxlen = 20000.0         # Upper limit for mixing length [m]
maxiso = 20000.0         # Upper limit for isotropy time scale [s]
w3clip = 1.2             # Third moment of vertical velocity clip
ustar_min = 0.01         # Minimum surface friction velocity
largeneg = -99999999.99  # Large negative value used for linear_interp threshold
dothetal_skew = False    # Temperature skewness independent of moisture variance
pblmaxp = 4e4            # PBL max depth in pressure units [Pa]
