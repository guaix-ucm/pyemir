import matplotlib

pytest_plugins = ["numina.testing.plugins"]

# The recipes draw plots, use a backend without display
matplotlib.use("Agg")
