from pathlib import Path
from typing import Tuple
import re

def version() -> Tuple[int]:
    version_file = Path(__file__).parent.parent.resolve() / "VERSION"

    if not version_file.exists():
        raise FileNotFoundError("Missing VERSION file in MMA root")

    return tuple(
        map(
            lambda x: int(x),
            re.match(
                r"^(\d+).(\d+).(\d+)",
                version_file.read_text()
            ).groups()
        )
    )

_version = version()