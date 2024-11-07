"""Console script for gdb_brexazine."""
import gdb_brexazine

import typer
from rich.console import Console

app = typer.Typer()
console = Console()


@app.command()
def main():
    """Console script for gdb_brexazine."""
    console.print("Replace this message by putting your code into "
               "gdb_brexazine.cli.main")
    console.print("See Typer documentation at https://typer.tiangolo.com/")
    


if __name__ == "__main__":
    app()
