"""JAX ports of E3SMv3 (EAM) parameterizations. float64 everywhere."""
import jax

jax.config.update("jax_enable_x64", True)
