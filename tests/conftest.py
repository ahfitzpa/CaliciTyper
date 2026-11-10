# Import calicityper before anything pulls in NumPy, so its BLAS thread pin
# (calicityper/__init__.py) is in place before NumPy starts a thread pool.
# Without this, a test module that imports Biopython first gets one BLAS
# thread per core, which can exhaust per-user thread limits on shared nodes.
import calicityper  # noqa: F401
