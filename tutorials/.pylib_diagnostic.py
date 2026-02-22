modules = [
"matplotlib", 
"scipy", 
"numpy", 
"os",
"sys",
"shutil",
"subprocess",
"time"]

missing = []

for module in modules:
    try:
        __import__(module)
    except ImportError:
        missing.append(module)

if missing:
    raise ImportError(
        "Missing required libraries: " + ", ".join(missing)
    )

print("all required python libraries are installed")

import sys

min_version = (3, 2)

if sys.version_info < min_version:
    raise RuntimeError(
        f"python {min_version[0]}.{min_version[1]} or newer is required. "
        f"you are using {sys.version.split()[0]}"
    )
