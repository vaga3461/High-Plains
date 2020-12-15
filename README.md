# Welcome

This document will be updated regularly to include most current information as this project progresses. This is also where you can find examples (after the model actually works)

# 12/10/2020 Update
- Having troubles with numpy not being recognized within my model class (getting `'np' is not defined`)

> I'm guessing this was because originally you did not have numpy imported in the `model.py` file and then it did not auto-reload without kernel restart (see note on enabling auto reload to address this issue). See other changes for how to fix the class.

- When trying to use function `run_model`, am getting `'self' is not defined`

> I cannot figure out what `run_model` is supposed to do. It seems it is initalizing another model inside the model object, which is probably not what you want to do. Added a `run_model()` example to show how to call it in CT_1D_erosion.ipynb, take a look.

# Goals
- To write a model that calculates erosion and deposition across a 1D river profile
- To encapsulate that model as a Python class within a function that can be imported into a notebook (this means that a user will not have to see all the messy, under-the-hood code)
- To define functions that represent different styles of tectonic uplift and tilting, which can also be imported into a notebook

## How to Use
- Import an uplift function from `functions.py`
- Import the model from `model.py`
- Only requires basic Python libraries (numpy and matplotlib)
