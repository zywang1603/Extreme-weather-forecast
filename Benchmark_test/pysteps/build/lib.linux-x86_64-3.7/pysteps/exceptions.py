# -*- coding: utf-8 -*-

# Custom pySteps exceptions


class MissingOptionalDependency(Exception):
    """Raised when an optional dependency is needed but not found."""
    pass


class DataModelError(Exception):
    """Raised when a file is not compilant with the Data Information Model."""
    pass
