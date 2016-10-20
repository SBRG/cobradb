class AlreadyLoadedError(Exception):
    pass

from cobradb.loading.model_loading import load_model
from cobradb.loading.component_loading import load_genome
