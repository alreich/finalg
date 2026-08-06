
import json

from finalg.make_finite_algebra import make_finite_algebra

class Examples:
    """A convenience class for retrieving some example algebras in the "algebras"
    directory.  To add or subtract algebras to its default list, see the file,
    'examples.json', in the "algebras" directory.

    ``algebras_dir`` may be a plain filesystem path (str/os.PathLike) or an
    ``importlib.resources`` Traversable (as returned by
    ``importlib.resources.files(...)``) -- the latter is what finalg's
    ``__init__.py`` passes in, so that example algebras are located
    correctly whether finalg is run from a source checkout or installed
    from a wheel."""

    def __init__(self, algebras_dir, filenames_json='examples.json'):
        examples_path = algebras_dir / filenames_json
        with examples_path.open('r') as fin:
            self.filenames_list = json.load(fin)
        self.algebras = [make_finite_algebra(str(algebras_dir / filename))
                         for filename in self.filenames_list]
        # self.about()

    def __len__(self):
        return len(self.algebras)

    def __getitem__(self, index):
        return self.algebras[index]

    def about(self):
        """Returns a list of example algebras with instructions on how to retrieve them."""
        n = 70
        print("=" * n)
        print(" " * (int(n / 2) - 8) + "Example Algebras")  # centered text
        print("-" * n)
        print(f"  {len(self.algebras)} example algebras are available.")
        # print("  Use \"Examples[INDEX]\" to retrieve a specific example,")
        print("  The INDEX is the first number on each line below:")
        print("-" * n)
        index = 0
        for alg in self.algebras:
            print(f"{index}: {alg.name} -- {alg.description}")
            index += 1
        print("=" * n)


