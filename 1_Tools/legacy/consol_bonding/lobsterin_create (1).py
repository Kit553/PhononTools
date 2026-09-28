"""Compatibility wrapper for the package-native lobsterin generator."""

try:
    from .lobsterin import main
except ImportError:  # Allows ``python path/to/lobsterin_create.py ...`` after installation.
    from cat_tools.bonding.lobsterin import main

if __name__ == "__main__":
    main()
