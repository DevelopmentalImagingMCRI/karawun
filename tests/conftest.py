import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--runcreatebaseline",
        action="store_true",
        default=False,
        help="create sha file"
    )


def pytest_configure(config):
    config.addinivalue_line("markers",
                            ("createbaseline: mark test as creating "
                             "baseline data"))


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runcreatebaseline"):
        # --runcreatebaseline given in cli: do not skip slow tests
        return
    skip_createbaseline = pytest.mark.skip(
        reason="need --runcreatebaseline option to run")
    for item in items:
        if "createbaseline" in item.keywords:
            item.add_marker(skip_createbaseline)
