import os
import json
import datetime


class InextensibleDict(dict):
    """A dictionary that refuses adding new keys -- safe against common typos."""
    def __setitem__(self, key, item):
        if key not in self:
            raise UserParameterError(
                "not possible to add new key {!r} to the parameters."
                .format(key)
            )
        dict.__setitem__(self, key, item)


class UserParameterError(Exception):
    """An error occured with the user-set parameters."""


class DefaultsFileError(Exception):
    """An error occured in the file of default values."""
    def __init__(self, defaults_file, arg):
        self.args = ["bad format of file {}: {}".format(defaults_file, arg)]


class UserParameters(object):
    """User interface for modifing the experiment parameters.

    The default parameters are stored in a JSON file in the core-
    directory.  Its filename is stored in the constant attribute
    DEFAULTS_FILE of this class.

    Attributes for public access:
     - model, IO, time, discretization, MPI: dictionaries with the
        parameters of the respective category.  It is possible to modify
        the value of each parameter, but it is not possible to add new
        parameters to any of the dictionaries.

    Methods for public access:
     - help
     - possible_values
     - view_parameters
     - check
    """

    DEFAULTS_FILE = "defaults.json"

    # Dictionary of datatypes; used to check that every parameter is of
    # correct type; key is a string as the type appears in the file of
    # default values; value is a Python type or a list of Python types
    TYPES = {
        "str": str,
        "int": int,
        "float": (float, int),  # an int can be used instead of a float
        "bool": bool,
        "list or string": (list, tuple, str),
    }

    def __init__(self):
        """Load the parameters from the file of default values."""
        # Get the full path of the file of default values
        jsonfile = os.path.realpath(os.path.join(
            os.getcwd(),
            os.path.dirname(__file__),
            self.DEFAULTS_FILE,
        ))
        # Read and check the file of default values
        with open(jsonfile) as f:
            defaults = json.load(f)
        self.check_defaults(defaults)

        # Copy all parameters by category
        self.model = InextensibleDict({
            parameter: attributes["default"]
            for parameter, attributes in defaults["model"].items()
        })
        self.IO = InextensibleDict({
            parameter: attributes["default"]
            for parameter, attributes in defaults["IO"].items()
        })
        self.animation = InextensibleDict({
            parameter: attributes["default"]
            for parameter, attributes in defaults["animation"].items()
        })
        self.physics = InextensibleDict({
            parameter: attributes["default"]
            for parameter, attributes in defaults["physics"].items()
        })
        self.time = InextensibleDict({
            parameter: attributes["default"]
            for parameter, attributes in defaults["time"].items()
        })
        self.discretization = InextensibleDict({
            parameter: attributes["default"]
            for parameter, attributes in defaults["discretization"].items()
        })
        self.MPI = InextensibleDict({
            parameter: attributes["default"]
            for parameter, attributes in defaults["MPI"].items()
        })

        # Copy information used for help and checks
        self.documentations = {}
        self.options = {}
        self.types = {}
        for category, parameters in defaults.items():
            for parameter, attributes in parameters.items():
                self.documentations[parameter] = attributes["doc"]
                self.options[parameter] = attributes["avail"]
                self.types[parameter] = attributes["type"]

    def help(self, parameter):
        """Return the documentation for a parameter."""
        if parameter in self.documentations:
            return self.documentations[parameter]
        raise ValueError("invalid parameter: {!r}".format(parameter))

    def possible_values(self, parameter):
        """Return the possible values for a parameter (list or string)."""
        if parameter in self.options:
            return self.options[parameter]
        raise ValueError("invalid parameter: {!r}".format(parameter))

    def view_parameters(self):
        """Have a look at all of the parameters in one dictionary.

        Warning: it is not possible to modify parameters using the
        dictionary returned by this function.
        """
        return {
            **self.model,
            **self.physics,
            **self.IO,
            **self.animation,
            **self.time,
            **self.discretization,
            **self.MPI,
        }

    def check(self):
        """Raise an exception if the value of any parameter is wrong."""
        # List of values that are recognized as powers of two.  Rationale:
        # In 2019, it seems ridicoulous to use 2^19 for any parameter in the
        # model, so this number raises a "too large" error.  Assuming that
        # this limit doubles every year (which is faster growth than
        # predicted by Moore's law), the following list, which adds one more
        # power of two every year, should be safe for the near future.
        POWERS_OF_2 = [2**n for n in range(datetime.datetime.now().year - 2000)]
        for parameter, value in self.view_parameters().items():
            # Check type of value
            param_type = self.types[parameter]
            if not isinstance(value, self.TYPES[param_type]):
                raise UserParameterError(
                    "parameter {} must be a {}, not {}"
                    .format(parameter, param_type, type(value))
                )
            # Check if value is among the options
            options = self.options[parameter]
            if isinstance(options, list):
                if not value in options:
                    raise UserParameterError(
                        "parameter {} must be one of {}, not {!r}"
                        .format(parameter, options, value)
                    )
            elif options == "> 0.0":
                if not value > 0.0:
                    raise UserParameterError(
                        "parameter {} must be positive".format(parameter)
                    )
            elif options == ">= 0.0" or options == ">= 0":
                if not value >= 0.0:
                    raise UserParameterError(
                        "parameter {} must be non-negative".format(parameter)
                    )
            elif options == "2^n":
                if value not in POWERS_OF_2:
                    if value < max(POWERS_OF_2):
                        raise UserParameterError(
                            "parameter {} must be a power of 2"
                            .format(parameter)
                        )
                    else:
                        raise UserParameterError(
                            "parameter {} is very large; if you are sure to "
                            "use it, extend the variable POWERS_OF_2 in this "
                            "class to include your value".format(parameter)
                        )
            elif options == "any valid filename":
                if "/" in value:
                    raise UserParameterError(
                        'parameter {} must not contain a "/"'.format(parameter)
                    )
            elif parameter == "variables_in_history":
                if not isinstance(value, (list, tuple)) and value != "all":
                    raise UserParameterError(
                        'parameter {} must be a list or "all", not {!r}'
                        .format(parameter, value)
                    )
            elif parameter == "datadir":
                # no check necessary
                pass
            else:
                # this should not happen; if it does, modify this method
                print("Warning: cannot check parameter", parameter)
        # Check that every CPU has at least one grid point
        for x in "xyz":
            if self.MPI["np" + x] > self.discretization["global_n" + x]:
                raise UserParameterError(
                    "parameter np{} cannot be larger than global_n{}"
                    .format(x, x)
                )

    @classmethod
    def check_defaults(cls, defaults):
        """Raise an exception if an error with the defaults is found.

        This does not check if the given default value is valid."""
        # List of categories that are needed in the defaults
        CATEGORIES = ["model", "physics", "IO", "animation", "time", "discretization", "MPI"]
        # List of attributes that are needed for every parameter
        ATTRIBUTES = ["type", "default", "avail", "doc"]
        # Check all categories exist
        for category in CATEGORIES:
            if category not in defaults:
                raise DefaultsFileError(
                    cls.DEFAULTS_FILE,
                    "category {} missing".format(category),
                )
        parameters = []
        for category in defaults:
            # Check no additional category exist
            if category not in CATEGORIES:
                raise DefaultsFileError(
                    cls.DEFAULTS_FILE,
                    "unknown category {}".format(category),
                )
            for parameter, attributes in defaults[category].items():
                # Check no parameter exists twice
                if parameter in parameters:
                    raise DefaultsFileError(
                        cls.DEFAULTS_FILE,
                        "parameter {} appears twice".format(parameter),
                    )
                parameters.append(parameter)
                # Check every parameter has all attributes
                for attribute in ATTRIBUTES:
                    if attribute not in attributes:
                        raise DefaultsFileError(
                            cls.DEFAULTS_FILE,
                            "attribute {} missing for parameter {}"
                            .format(attribute, parameter),
                        )
                # Check no parameter has additional attributes
                for attribute in attributes:
                    if attribute not in ATTRIBUTES:
                        raise DefaultsFileError(
                            cls.DEFAULTS_FILE,
                            "unknown attribute {} for parameter {}"
                            .format(attribute, parameter),
                        )
                # Check parameter type is recognized
                if attributes["type"] not in cls.TYPES:
                    raise DefaultsFileError(
                        cls.DEFAULTS_FILE,
                        "unknown type {} of parameter {}"
                        .format(attributes["type"], parameter),
                    )


if __name__ == "__main__":
    param = UserParameters()

    print("Default parameters:", param.view_parameters())
    print("-"*80)

    # Changing a parameter
    param.IO["expname"] = "2nd_experiment"
    # Changing a parameter that does not exist
    try: param.IO["epxname"] = "3rd_experiment"  # intentional typo
    except UserParameterError as e: print("UserParameterError:", e)
    # Changing a parameter in the wrong set
    try: param.model["expname"] = "4th_experiment"
    except UserParameterError as e: print("UserParameterError:", e)
    print("-"*80)

    print("Model parameters:", param.model)
    print("Physics parameters:", param.physics)
    print("Input/Output parameters:", param.IO)
    print("Animation parameters:", param.animation)
    print("Time parameters:", param.time)
    print("Discretization parameters:", param.discretization)
    print("MPI parameters:", param.MPI)
    print("-"*80)

    print("Possible values for modelname:", param.possible_values("modelname"))
    print("Description of modelname:", param.help("modelname"))
    print("-"*80)

    param.check()

    # Count the exceptions to make it easy to check if all were raised
    i = 0

    # Some examples that create a failing check
    i += 1
    param.model["geometry"] = "open"
    try: param.check()
    except UserParameterError as e: print("{:2d}. UserParameterError: {}".format(i, e))
    param.model["geometry"] = "closed"

    i += 1
    param.model["Lx"] = -3
    try: param.check()
    except UserParameterError as e: print("{:2d}. UserParameterError: {}".format(i, e))
    param.model["Lx"] = 3

    i += 1
    param.model["Ly"] = "5"
    try: param.check()
    except UserParameterError as e: print("{:2d}. UserParameterError: {}".format(i, e))
    param.model["Ly"] = 5

    i += 1
    param.IO["expname"] = "foo/bar"
    try: param.check()
    except UserParameterError as e: print("{:2d}. UserParameterError: {}".format(i, e))
    param.IO["expname"] = "foo-bar"

    i += 1
    param.IO["variables_in_history"] = "any"
    try: param.check()
    except UserParameterError as e: print("{:2d}. UserParameterError: {}".format(i, e))
    param.IO["variables_in_history"] = "all"

    i += 1
    param.IO["timestep_history"] = -1.0
    try: param.check()
    except UserParameterError as e: print("{:2d}. UserParameterError: {}".format(i, e))
    param.IO["timestep_history"] = 0.0

    i += 1
    param.time["cfl"] = 0
    try: param.check()
    except UserParameterError as e: print("{:2d}. UserParameterError: {}".format(i, e))
    param.time["cfl"] = 1

    i += 1
    param.discretization["global_nx"] = 3
    try: param.check()
    except UserParameterError as e: print("{:2d}. UserParameterError: {}".format(i, e))
    param.discretization["global_nx"] = 2

    i += 1
    param.discretization["global_ny"] = 2**1000
    try: param.check()
    except UserParameterError as e: print("{:2d}. UserParameterError: {}".format(i, e))
    param.discretization["global_ny"] = 2**10

    i += 1
    param.discretization["global_nz"] = 1
    param.MPI["npz"] = 2
    try: param.check()
    except UserParameterError as e: print("{:2d}. UserParameterError: {}".format(i, e))
    param.MPI["npz"] = 1

    i += 1
    param.MPI["npx"] = 1.5
    try: param.check()
    except UserParameterError as e: print("{:2d}. UserParameterError: {}".format(i, e))
    param.MPI["npx"] = 1

    i += 1
    param.MPI["nh"] = -1
    try: param.check()
    except UserParameterError as e: print("{:2d}. UserParameterError: {}".format(i, e))
    param.MPI["nh"] = 1

    i += 1
    param.discretization["orderA"] = 4
    try: param.check()
    except UserParameterError as e: print("{:2d}. UserParameterError: {}".format(i, e))
    param.discretization["orderA"] = 3
