class PythonOwnedError(Exception):
    """Error class to use in case it is forbidden to use Python classes"""
    def __init__(
        self,
        message="Python initialized ctypes structure cannot be used as argument. "
    ):
        super().__init__(message)

class InitializedStructureError(Exception):
    """Error class to use in case it is forbidden to use initialized ctypes classes"""
    def __init__(
        self,
        message=(
            "Initialized ctypes structure with allocated fields "
            "cannot be used as argument."
        )
    ):
        super().__init__(message)

class NotInitializedError(RuntimeError):
    """Raised if DLL is not set up before ctypes structure instantiation"""
    def __init__(
        self,
        message=(
            "Python TDSE DLL has not been initialized yet. "
            "Create an instance of TDSE_DLL class first."
        )
    ):
        super().__init__(message)