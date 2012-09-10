

EMIR Data Products
==================

Each recipe of the EMIR Pipeline produces a set of predefined results, known
as *data products*. In turn, the recipes may request diferent data products
as computing requeriments, efectively chaining the recipes.

For example, the requirements of the :ref:`ff-recipe-label` recipe include a
:class:`~emir.dataproducts.MasterDark` object. This object is produced by 
the recipe :class:`~emir.recipes.DarkRecipe`, which in turn requires a 
:class:`~emir.dataproducts.MasterBias` object.


