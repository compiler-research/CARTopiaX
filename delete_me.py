import numpy
import pandas

try:
    import optuna
    optuna_version = optuna.__version__
except ImportError:
    optuna_version = "NO INSTALADO"

print("numpy:", numpy.__version__)
print("pandas:", pandas.__version__)
print("optuna:", optuna_version)