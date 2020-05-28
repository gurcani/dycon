# DyCoN Models

Dynamical Complex Network models as described in the paper https://arxiv.org/abs/2004.01134.

Each run_XXX directory contains a dycon_gns.py, which is basically the same file but with different parameters so that it corresponds to a different variant of the model. You should run all the models (which may take up to a day in total) before you can run the plot routines such as plot_eval.py etc.

On the other hand routines such as draw_bipartite_nw.py or plot_dd.py depends only on static networks found in corresponding model directories, thus one does not need to run the models to see features of the static networks.
