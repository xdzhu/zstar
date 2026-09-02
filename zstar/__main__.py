"""Allow ``python -m zstar`` to use the installed command-line interface."""

from .cli import zstar_cli


if __name__ == "__main__":
    zstar_cli()
